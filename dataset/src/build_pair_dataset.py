"""Build the PAIR longitudinal master table and the v0.1 dataset release.

Row definition: one patient x one longitudinal sample interval x one treatment context.
Everything is deterministically derived from the canonical corpus parquet tables plus
the explicit recovery sources in corpus/src/patient_archive.py. No second truth source
is created; this script only projects existing records into the release layout.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path

import pandas as pd

DATASET_ROOT = Path(__file__).resolve().parents[1]
CORPUS = DATASET_ROOT.parent / 'corpus'
PROJECT = DATASET_ROOT.parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import normalization as nz  # noqa: E402

MISSING = {'', 'NOT_FOUND', 'NOT_MEASURED', 'NOT_APPLICABLE', 'UNKNOWN', 'UNK', 'NA', 'nan'}
RESPONSE_FAMILIES = ['RECIST_RESPONSE', 'PCR', 'CHEMOTHERAPY_RESPONSE_SCORE',
                     'ENDPOINT_SEMANTICS_UNCERTAIN', 'LESION_LEVEL_RECIST_RESPONSE',
                     'CLINICAL_RESPONSE_ULTRASOUND_VOLUME', 'RECIST_BINARY_RESPONSE',
                     'CLINICAL_RESPONSE_COMPOSITE_NEOADJUVANT']
RECOVERED_FAMILIES = {'CHEMOTHERAPY_RESPONSE_SCORE', 'LESION_LEVEL_RECIST_RESPONSE',
                      'RECIST_BINARY_RESPONSE'}
# Same study, shared GSM sample accessions; resource-level duplicate whose patients
# are merged at patient_uid level on the GSM-confirmed crosswalk (see below).
SAME_STUDY_RESOURCE_GROUPS = {'DUP_GSE20181_GSE5462': ['GSE20181', 'GSE5462']}


def load_archive():
    spec = importlib.util.spec_from_file_location('patient_archive', CORPUS / 'src/patient_archive.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def portable_locator(value):
    """Turn a corpus file reference into a locator that is portable and still identifying.

    A published record may say where its source lives *inside this project*; it must not say
    where one person's machine keeps it.  Absolute paths inside the project tree become
    project-relative POSIX, which is the same convention ``build_gene_space.SOURCES`` uses, so
    a reader can resolve them with no configuration.  An already-relative value is taken as
    project-relative and normalised, never re-resolved against the working directory — doing
    that invents a path that changes with the directory the build happened to run from.
    Anything else keeps its provenance in the accession, ``source_table``, ``source_column``,
    ``source_url`` and status columns and is labelled rather than being rewritten into an
    invented path.
    """
    text = '' if value is None else str(value)
    if not text or text in MISSING:
        return text
    parts = []
    for piece in text.split('|'):
        piece = piece.strip()
        if not piece:
            continue
        path = Path(piece)
        if not path.is_absolute():
            parts.append(path.as_posix())
            continue
        try:
            parts.append(path.resolve().relative_to(PROJECT).as_posix())
        except (ValueError, OSError):
            parts.append('SOURCE_EXTERNAL')
    return '|'.join(parts) if parts else text


def expression_reference(value):
    """MASTER's expression ref: a portable locator, or NOT_AVAILABLE when there is no file.

    A corpus row with no expression file records a sentinel, not a path.  Running that
    sentinel through a path resolver resolves it against the current working directory and
    publishes an invented locator that changes with wherever the build happened to start, so
    the sentinel is decided here before any path arithmetic.
    """
    if not value or str(value) in MISSING:
        return 'NOT_AVAILABLE'
    return portable_locator(value)


def resolved_frames(archive):
    """Run the verified archive pipeline and keep the in-memory frames."""
    tables = archive.load_tables(CORPUS)
    tables = archive.repair_gide(tables, CORPUS)
    recovered = [archive.recover_crs(CORPUS), archive.recover_gse179994(CORPUS),
                 archive.recover_gse120575(CORPUS)]
    endpoints, resolved_counts = archive.resolve_native_response_semantics(tables['clinical_endpoints'])
    endpoints = pd.concat([endpoints, *recovered], ignore_index=True).fillna('')
    author = pd.read_csv(PROJECT / 'data/raw/MORRISON-1-public/RNASeq/RNA-CancerCell-MORRISON1-metadata.tsv',
                         sep='\t').fillna('')
    samples = archive.resolve_samples(tables['samples'], author)
    endpoints = archive.link_endpoints(endpoints, samples)
    endpoints['source_path_status'] = endpoints.source_file.map(
        lambda value: 'LOCAL_SOURCE_PRESENT' if value and all(Path(p).exists() for p in value.split('|'))
        else 'SOURCE_PATH_UNAVAILABLE')
    transitions = tables['transitions'].copy()
    sample_index = samples.set_index('sample_id', verify_integrity=True)
    transitions['patient_key'] = transitions.sample_t0.map(sample_index.patient_key)
    end_key = transitions.sample_t1.map(sample_index.patient_key)
    transitions['pair_valid'] = [isinstance(a, str) and a == b for a, b in zip(transitions.patient_key, end_key)]
    exposures = tables['treatment_exposures'].merge(
        samples[['sample_id', 'patient_key']], on='sample_id', how='inner')
    gide_author = author[author.cohort.eq('gide')].set_index('sample.id')
    gide_rows = []
    for row in samples[samples.cohort_code.eq('GIDE')].itertuples():
        sid = row.repository_sample_id + '_gide'
        if sid in gide_author.index:
            gide_rows.append(dict(exposure_id=row.sample_id + ':treatment', sample_id=row.sample_id,
                                   patient_key=row.patient_key,
                                   native_treatment=gide_author.loc[sid, 'treatment.regimen.name'],
                                   treatment_id='', treatment_classes='', treatment_arm=''))
    exposures = pd.concat([exposures, pd.DataFrame(gide_rows)], ignore_index=True).fillna('')
    return tables, samples, endpoints, transitions, exposures, resolved_counts


def patient_index(archive, tables, samples, endpoints, transitions):
    """One row per resource-level patient entry with identity and endpoint summary."""
    gpi = tables['global_patient_identity']
    uid_lookup = dict(zip(gpi.local_key, gpi.union_root))
    conf_lookup = dict(zip(gpi.local_key, gpi.confidence))
    study_group = {cohort: gid for gid, cohorts in SAME_STUDY_RESOURCE_GROUPS.items() for cohort in cohorts}
    prior_lookup = {f'{r.cohort_code}:{r.native_patient_id}': r._asdict()
                    for r in tables['patients'].itertuples(index=False)}
    rows = []
    for key, group in samples.groupby('patient_key', sort=True):
        cohort, native = key.split(':', 1)
        legacy_keys = [cohort + ':' + str(n) for n in group.legacy_native_patient_id.unique()]
        prior = next((prior_lookup[k] for k in legacy_keys if k in prior_lookup), {})
        if cohort in ('GIDE', 'MORRISON_gide'):
            # Author-corrected identity; the legacy gpi root 'GIDE:13' wrongly merges the
            # two author-separated subjects PD1_13 and ipiPD1_13, so it is never used here.
            uid, conf = str(group.identity_group.iloc[0]), 'AUTHOR_SUBJECT_ID'
        else:
            uid, conf = '', 'UNKNOWN'
            for candidate in [key] + legacy_keys:
                if candidate in uid_lookup:
                    uid, conf = uid_lookup[candidate], conf_lookup[candidate]
                    break
        if not uid:
            uid = key
        records = endpoints[endpoints.patient_key.eq(key)]
        response = records[records.endpoint_family.isin(RESPONSE_FAMILIES)]
        status = archive.label_status(response.native_value.tolist(), response.endpoint_family.tolist())
        paired = transitions[transitions.patient_key.eq(key) & transitions.pair_valid]
        rows.append(dict(
            patient_key=key, patient_uid=uid, cohort_code=cohort, native_patient_id=native,
            identity_status=';'.join(sorted(set(group.identity_status))),
            identity_confidence=conf,
            duplicate_resource_group=study_group.get(cohort, ''),
            cancer_type=prior.get('cancer_type', 'NOT_ESTABLISHED'),
            cancer_type_native=prior.get('cancer_type_native', 'NOT_ESTABLISHED'),
            oncotree=prior.get('oncotree', 'NOT_MEASURED'),
            n_samples=len(group), has_pair=bool(len(paired)),
            endpoint_status=status,
            endpoint_families='|'.join(sorted(set(response.endpoint_family))),
            endpoint_native_values='|'.join(sorted({str(v) for v in response.native_value if str(v) not in MISSING})),
            endpoint_harmonized=nz.harmonize_endpoint(response.native_value.tolist(),
                                                      response.endpoint_family.tolist()),
            endpoint_assessment_time='|'.join(sorted({str(v) for v in response.assessment_time
                                                      if str(v) not in MISSING})),
            n_endpoint_rows=len(records),
            endpoint_verification='|'.join(
                f'{s}:{n}' for s, n in records.join_status.value_counts().items()) or 'NOT_LINKED',
            additive_universe=bool(prior.get('additive_universe', False))))
    patients = pd.DataFrame(rows)
    # Cross-cohort confirmed same-patient groups (GIDE/MORRISON_gide, GSE91061/MORRISON_038).
    shared = patients.patient_uid.duplicated(keep=False) & ~patients.patient_uid.eq('')
    for uid in patients.loc[shared, 'patient_uid'].unique():
        mask = patients.patient_uid.eq(uid) & patients.duplicate_resource_group.eq('')
        patients.loc[mask, 'duplicate_resource_group'] = uid
    # GSE20181/GSE5462 same-study resources: patients linked by shared GSM sample
    # accessions (116 shared, native ids agree on all) are the same biological
    # patients — merge at patient_uid level to prevent cross-resource train/test
    # leakage, while both accession-specific resource records are retained.
    crosswalk = gse20181_gse5462_crosswalk(samples)
    overlap = crosswalk[crosswalk.evidence.eq('SHARED_GSM_SAMPLE_ACCESSION')]
    merged_pairs = 0
    for left, right in overlap[['patient_gse20181', 'patient_gse5462']].drop_duplicates().itertuples(index=False):
        lkey, rkey = f'GSE20181:{left}', f'GSE5462:{right}'
        if lkey in patients.patient_key.values and rkey in patients.patient_key.values:
            patients.loc[patients.patient_key.eq(rkey), 'patient_uid'] = patients.loc[
                patients.patient_key.eq(lkey), 'patient_uid'].iloc[0]
            merged_pairs += 1
    patients.attrs['gsm_confirmed_merged_pairs'] = merged_pairs
    return patients


def build_master(archive, patients, samples, endpoints, transitions, exposures, treatments):
    treatment_classes = dict(zip(treatments.treatment_id, treatments.classes))
    native_classes = dict(zip(treatments.native, treatments.classes))
    sample_index = samples.set_index('sample_id', verify_integrity=True)
    exposure_by_sample = {}
    for sid, group in exposures.groupby('sample_id'):
        values = [v for v in group.native_treatment if v not in ('', 'NOT_MEASURED')]
        exposure_by_sample[sid] = '|'.join(sorted(set(values)))
    classes_by_sample = {}
    for sid, group in exposures.groupby('sample_id'):
        classes = ''
        for exp in group.itertuples(index=False):
            if exp.treatment_id in treatment_classes:
                classes = treatment_classes[exp.treatment_id]
                break
            if exp.native_treatment in native_classes:
                classes = native_classes[exp.native_treatment]
                break
        classes_by_sample[sid] = classes
    # Corpus sample ids addressable by each endpoint's native sample identifier.
    sample_ids_by_native = {}
    for row in samples.itertuples():
        for value in (row.native_sample_id, row.repository_sample_id):
            if value and value not in MISSING:
                sample_ids_by_native.setdefault((row.cohort_code, value), set()).add(row.sample_id)
    response_by_patient = {key: group for key, group in
                           endpoints[endpoints.endpoint_family.isin(RESPONSE_FAMILIES)]
                           .groupby('patient_key')}
    n_intervals_by_patient = transitions[transitions.pair_valid].groupby('patient_key').size()

    def interval_endpoints(key, sample_t1, cohort):
        """Response endpoints attributable to this interval, or (None, reason).

        Binding rules (fail-closed): an endpoint binds when it was measured on the
        t1 sample of this interval, or when the patient has exactly one valid
        interval so the attribution is unambiguous. Anything else stays at patient
        level only and is reported as NOT_INTERVAL_RESOLVED.
        """
        records = response_by_patient.get(key)
        if records is None or not len(records):
            return None, 'LABEL_NOT_FOUND'
        on_t1 = records[[sample_t1 in sample_ids_by_native.get((cohort, sid), set())
                         for sid in records.sample_id_native]]
        if len(on_t1):
            return on_t1, 'BOUND'
        if n_intervals_by_patient.get(key, 0) == 1:
            return records, 'BOUND'
        return None, 'NOT_INTERVAL_RESOLVED'

    patient_by_key = patients.set_index('patient_key')
    rows = []
    for tr in transitions.itertuples(index=False):
        key = tr.patient_key if isinstance(tr.patient_key, str) else ''
        patient = patient_by_key.loc[key] if key in patient_by_key.index else None
        s0, s1 = sample_index.loc[tr.sample_t0], sample_index.loc[tr.sample_t1]
        treatment_at_t1 = exposure_by_sample.get(tr.sample_t1, '')
        drug, drug_class, family = nz.classify_treatment(treatment_at_t1,
                                                         classes_by_sample.get(tr.sample_t1, ''))
        phase_t0, phase_t1 = nz.normalize_phase(s0.timepoint_normalized), nz.normalize_phase(s1.timepoint_normalized)
        # Treatment recorded "before the t1 sample" is attributable to this interval only
        # when it started after a PRE_TREATMENT t0; otherwise chronology is not established.
        if phase_t0 == 'PRE_TREATMENT' and phase_t1 in ('ON_TREATMENT_EARLY', 'ON_TREATMENT_LATE'):
            interval_treatment, treatment_status = treatment_at_t1, 'INTERVAL_ATTRIBUTED'
        else:
            interval_treatment, treatment_status = '', 'INTERVAL_ATTRIBUTION_UNRESOLVED'
        if patient is not None:
            p_status = patient.endpoint_status
            p_families = [f for f in str(patient.endpoint_families).split('|') if f]
            p_values = [v for v in str(patient.endpoint_native_values).split('|') if v]
            p_harmonized = str(patient.endpoint_harmonized)
            timing = str(patient.endpoint_assessment_time)
            verification = patient.endpoint_verification
            uid, identity, dupe = patient.patient_uid, patient.identity_status, patient.duplicate_resource_group
            cancer, cancer_native, subtype = patient.cancer_type, patient.cancer_type_native, patient.oncotree
            bound, reason = interval_endpoints(key, tr.sample_t1, tr.cohort_code)
            if reason == 'BOUND':
                i_status = archive.label_status(bound.native_value.tolist(), bound.endpoint_family.tolist())
                i_values = sorted({str(v) for v in bound.native_value if str(v) not in MISSING})
                i_families = sorted(set(bound.endpoint_family))
                i_timing = sorted({str(v) for v in bound.assessment_time if str(v) not in MISSING})
                i_harmonized = nz.harmonize_endpoint(bound.native_value.tolist(),
                                                     bound.endpoint_family.tolist())
            else:
                i_status, i_values, i_families, i_timing = reason, [], [], []
                i_harmonized = ''
        else:
            p_status, p_families, p_values, p_harmonized = 'MAPPING_NOT_FOUND', [], '', ''
            timing, verification = '', 'NOT_LINKED'
            uid, identity, dupe = key or 'UNMAPPED', 'MAPPING_NOT_FOUND', ''
            cancer = cancer_native = 'NOT_ESTABLISHED'
            subtype = 'NOT_MEASURED'
            i_status, i_values, i_families, i_timing = 'MAPPING_NOT_FOUND', [], [], []
        rows.append(dict(
            pair_uid=tr.transition_id,
            cohort_code=tr.cohort_code,
            patient_uid=uid,
            native_patient_id=key.split(':', 1)[1] if ':' in key else '',
            identity_status=identity,
            duplicate_resource_group=dupe,
            cancer_type_native=cancer_native,
            cancer_type=cancer,
            disease_subtype=subtype,
            sample_t0=tr.sample_t0, sample_t1=tr.sample_t1,
            time_t0_native=tr.time_t0_native, time_t1_native=tr.time_t1_native,
            phase_t0=phase_t0, phase_t1=phase_t1,
            delta_days='NOT_MEASURED',
            interval_treatment=interval_treatment or 'NOT_MEASURED',
            interval_treatment_status=treatment_status,
            treatment_at_t1=treatment_at_t1 or 'NOT_MEASURED',
            drug=drug, drug_class=drug_class, treatment_family=family,
            patient_endpoint_family_native='|'.join(p_families),
            patient_endpoint_family_harmonized=nz.endpoint_families_master(p_families),
            patient_endpoint_system=nz.endpoint_systems(p_families),
            patient_endpoint_native='|'.join(p_values),
            patient_endpoint_harmonized=p_harmonized,
            patient_endpoint_assessment_time=timing,
            patient_endpoint_status=p_status,
            patient_endpoint_verification_status=verification,
            interval_endpoint_status=i_status,
            interval_endpoint_family_native='|'.join(i_families),
            interval_endpoint_family_harmonized=nz.endpoint_families_master(i_families),
            interval_endpoint_system=nz.endpoint_systems(i_families),
            interval_endpoint_native='|'.join(i_values),
            interval_endpoint_harmonized=i_harmonized,
            interval_endpoint_assessment_time='|'.join(i_timing),
            expression_modality=s1.modality or 'NOT_MEASURED',
            platform=s1.platform or 'NOT_MEASURED',
            expression_t0_ref=expression_reference(s0.expression_file),
            expression_t1_ref=expression_reference(s1.expression_file),
            gene_space=nz.gene_space_for(s1.modality),
            source_accession=s1.dataset_id or tr.cohort_code,
            source_publication=nz.SOURCE_PUBLICATION.get(tr.cohort_code, ''),
            source_table='longitudinal-data/corpus/data/transitions.parquet',
            provenance_id=tr.transition_id,
            paired_valid=bool(tr.pair_valid),
            clinical_endpoint_available=p_status.startswith('OBSERVED_'),
            strict_prcr_vs_pd_eligible=bool(tr.pair_valid) and i_status == 'OBSERVED_STRICT_RECIST',
            frozen_auo_eligible=bool(tr.pair_valid) and i_status == 'OBSERVED_STRICT_RECIST'
            and tr.cohort_code == 'GSE91061',
        ))
    return pd.DataFrame(rows)


def gse20181_gse5462_crosswalk(samples):
    """Explicit sample-level crosswalk for the same-study duplicate resources.

    Links are established only by shared GSM sample accessions; patient ids are
    compared for information but never used to merge records.
    """
    columns = ['gsm', 'sample_gse20181', 'patient_gse20181', 'sample_gse5462', 'patient_gse5462',
               'native_patient_ids_agree', 'evidence']
    rows = []
    left = samples[samples.cohort_code.eq('GSE20181')].set_index('repository_sample_id', drop=False)
    right = samples[samples.cohort_code.eq('GSE5462')].set_index('repository_sample_id', drop=False)
    seen = set()
    for gsm in sorted(set(left.index) | set(right.index)):
        l = left.loc[gsm] if gsm in left.index else None
        r = right.loc[gsm] if gsm in right.index else None
        if l is not None and r is not None:
            agree = str(l.native_patient_id) == str(r.native_patient_id)
            rows.append(dict(zip(columns, [gsm, l.sample_id, l.native_patient_id, r.sample_id,
                                           r.native_patient_id, agree,
                                           'SHARED_GSM_SAMPLE_ACCESSION'])))
            seen.add(gsm)
        elif l is not None:
            rows.append(dict(zip(columns, [gsm, l.sample_id, l.native_patient_id, '', '',
                                           'NOT_COMPARABLE', 'GSM_PRESENT_ONLY_IN_GSE20181'])))
        else:
            rows.append(dict(zip(columns, [gsm, '', '', r.sample_id, r.native_patient_id,
                                           'NOT_COMPARABLE', 'GSM_PRESENT_ONLY_IN_GSE5462'])))
    return pd.DataFrame(rows, columns=columns)


def duplicate_groups(patients):
    rows = []
    for gid, group in patients.groupby('duplicate_resource_group', sort=True):
        if not gid:
            continue
        members = sorted(group.patient_key)
        uid_sizes = group.groupby('patient_uid').size()
        merged_uids = int((uid_sizes > 1).sum())
        merged_members = int(uid_sizes[uid_sizes > 1].sum())
        if len(uid_sizes) == 1:
            group_type = 'SAME_PATIENT'
            evidence = 'GLOBAL_PATIENT_IDENTITY_CONFIRMED_SAME_PATIENT'
        elif merged_uids and gid in SAME_STUDY_RESOURCE_GROUPS:
            group_type = 'SAME_STUDY_RESOURCE'
            evidence = 'SHARED_GSM_SAMPLE_ACCESSIONS_CONFIRMED_SAME_PATIENT'
        else:
            group_type = 'SAME_STUDY_RESOURCE'
            evidence = 'SHARED_GSM_SAMPLE_ACCESSIONS_SAME_STUDY'
        rows.append(dict(
            duplicate_group_id=gid,
            group_type=group_type,
            n_members=len(members), members='|'.join(members),
            evidence=evidence,
            merged_patient_uids=merged_uids,
            unmerged_members=len(group) - merged_members,
            patients_merged='YES' if merged_members == len(group) and merged_uids
            else ('PARTIAL' if merged_uids else 'NO')))
    return pd.DataFrame(rows)


def provenance_table(master, endpoints):
    endpoint_sources = {}
    for key, group in endpoints[endpoints.patient_key.ne('')].groupby('patient_key'):
        endpoint_sources[key] = '|'.join(sorted(
            {Path(p).name for v in group.source_file for p in str(v).split('|') if p}))
    rows = []
    for row in master.itertuples(index=False):
        key = f'{row.cohort_code}:{row.native_patient_id}'
        families = set(row.patient_endpoint_family_native.split('|'))
        rows.append(dict(
            pair_uid=row.pair_uid, patient_key=key, patient_uid=row.patient_uid,
            transition_source='corpus/data/transitions.parquet',
            endpoint_source_files=endpoint_sources.get(key, ''),
            endpoint_source='corpus/clinical/clinical_endpoints.parquet + recovery sources (see LABEL_RECOVERY_LEDGER.csv)',
            recovery_labels_used='YES' if families & RECOVERED_FAMILIES else 'NO'))
    return pd.DataFrame(rows)


def summarize(master, patients, samples, crosswalk):
    endpoint_verified = patients.endpoint_status.str.startswith('OBSERVED_')
    strict = master.paired_valid & master.interval_endpoint_status.eq('OBSERVED_STRICT_RECIST')
    frozen = master.frozen_auo_eligible
    per_cohort = []
    for cohort, group in patients.groupby('cohort_code', sort=True):
        master_rows = master[master.cohort_code.eq(cohort)]
        per_cohort.append(dict(
            cohort_code=cohort, resource_patient_entries=len(group),
            patient_uids=group.patient_uid.nunique(),
            master_rows=len(master_rows),
            paired_patients=int(group.has_pair.sum()),
            endpoint_verified_patients=int(endpoint_verified[group.index].sum()),
            paired_endpoint_verified=int((group.has_pair & endpoint_verified[group.index]).sum()),
            strict_prcr_vs_pd_eligible_patients=int(master_rows[strict[master_rows.index]]
                                                    .groupby('patient_uid').ngroups),
            cancer_type=group.cancer_type.iloc[0]))
    dupes = patients[patients.duplicate_resource_group.ne('')]
    shared = crosswalk[crosswalk.evidence.eq('SHARED_GSM_SAMPLE_ACCESSION')]
    overlap_pairs = shared.groupby(['patient_gse20181', 'patient_gse5462']).ngroups
    merged_pairs = patients.attrs.get('gsm_confirmed_merged_pairs', 0)
    same_study = [dict(
        duplicate_group_id='DUP_GSE20181_GSE5462', cohorts='GSE20181|GSE5462',
        group_type='SAME_STUDY_RESOURCE',
        patients_merged='YES' if merged_pairs == overlap_pairs else 'PARTIAL',
        merged_patient_pairs=int(merged_pairs),
        evidence='SHARED_GSM_SAMPLE_ACCESSIONS_CONFIRMED_SAME_PATIENT',
        shared_gsm_samples=int(len(shared)),
        patient_overlap_by_shared_gsm=int(overlap_pairs),
        note='both accession-specific resource records retained; biological patients '
             'merged at patient_uid level to prevent cross-resource train/test leakage')]
    return dict(
        pair_dataset_version='v0.1.2',
        generated_from='longitudinal-data/corpus canonical parquet tables + patient_archive recovery pipeline',
        cohort_count=len(per_cohort),
        resource_patient_entries=len(patients),
        patient_uids_after_confirmed_same_patient_merges=int(patients.patient_uid.nunique()),
        gsm_confirmed_merged_patient_pairs=int(merged_pairs),
        same_study_duplicate_resources=same_study,
        sample_records=len(samples),
        master_rows=len(master),
        longitudinal_pairs=len(master),
        paired_patients=int(patients.has_pair.sum()),
        cancer_types=int(patients.cancer_type.nunique()),
        treatment_families=sorted(f for f in master.treatment_family.unique() if f not in ('', 'NOT_MEASURED')),
        assay_modalities=sorted(m for m in samples.modality.unique() if m),
        endpoint_families=sorted(f for f in patients.endpoint_families.str.split('|').explode().unique() if f),
        endpoint_verified_patients=int(endpoint_verified.sum()),
        paired_endpoint_verified_patients=int((patients.has_pair & endpoint_verified).sum()),
        strict_prcr_vs_pd_eligible_patients=int(master[strict].groupby('patient_uid').ngroups),
        frozen_auo_eligible_patients=int(master[frozen].groupby('patient_uid').ngroups),
        interval_endpoint_status_counts=master.interval_endpoint_status.value_counts().to_dict(),
        duplicate_resource_groups=int(dupes.duplicate_resource_group.nunique()),
        by_cohort=per_cohort,
        note='patient_uids_after_confirmed_same_patient_merges counts resource entries merged '
             'only on confirmed same-patient evidence (GIDE/MORRISON_gide author subject IDs, '
             'GSE91061/MORRISON_038 confirmed roots, GSE20181/GSE5462 shared-GSM patients); '
             'it is NOT a final unique biological-patient count for cohorts lacking identity '
             'confirmation across resources (see same_study_duplicate_resources)')


SNAPSHOT_INPUTS = [
    # (path relative to the melanoma-immune-escape project root, public origin)
    ('longitudinal-data/corpus/clinical/clinical_endpoints.parquet', 'derived from public GEO/series metadata (see LABEL_RECOVERY_LEDGER.csv)'),
    ('longitudinal-data/corpus/data/patients.parquet', 'derived from public cohort metadata'),
    ('longitudinal-data/corpus/data/samples.parquet', 'derived from public cohort metadata'),
    ('longitudinal-data/corpus/data/transitions.parquet', 'derived from public cohort metadata'),
    ('longitudinal-data/corpus/data/treatments.parquet', 'derived from public cohort metadata'),
    ('longitudinal-data/corpus/data/assays.parquet', 'derived from public cohort metadata'),
    ('longitudinal-data/corpus/data/treatment_exposures.parquet', 'derived from public cohort metadata'),
    ('longitudinal-data/corpus/data/cohorts.parquet', 'derived from public cohort metadata'),
    ('longitudinal-data/corpus/data/exclusions.parquet', 'derived curation table'),
    ('longitudinal-data/corpus/data/timepoints.parquet', 'derived from public cohort metadata'),
    ('longitudinal-data/corpus/data/trajectories.parquet', 'derived from public cohort metadata'),
    ('longitudinal-data/corpus/data/warnings.parquet', 'derived curation table'),
    ('longitudinal-data/corpus/data/clinical_endpoints.parquet', 'derived from public cohort metadata (clinical/ copy is the authority)'),
    ('longitudinal-data/corpus/metadata/global_patient_identity.parquet', 'derived from public cohort metadata'),
    ('longitudinal-data/corpus/metadata/eligibility.parquet', 'derived curation table'),
    ('longitudinal-data/corpus/metadata/patient_aliases.parquet', 'derived curation table'),
    ('longitudinal-data/corpus/metadata/leakage_groups.parquet', 'derived curation table'),
    ('data/raw/MORRISON-1-public/RNASeq/RNA-CancerCell-MORRISON1-metadata.tsv', 'https://github.com/livnatje/MORRISON-1-public (author metadata)'),
    ('longitudinal-data/data/raw/GSE165897/PMC8865800.xml', 'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8865800/'),
    ('longitudinal-data/data/raw/GSE179994/supplement/43018_2021_292_MOESM3_ESM.xlsx', 'https://doi.org/10.1038/s43018-021-00292-8 (Supplementary Table 1)'),
    ('longitudinal-data/data/raw/GSE120575/GSE120575_series_matrix.txt.gz', 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120575'),
]


# The only two columns in this project that hold a file reference.  A redistributed snapshot
# must not carry a private machine layout, so these are converted on export; a future path
# column is caught by the published-artifact invariant test rather than by guessing here.
PORTABLE_LOCATOR_COLUMNS = ('expression_file', 'source_file')


def portable_frame(path: Path):
    """Read a redistributable table and rewrite its file locators in portable form."""
    frame = pd.read_parquet(path)
    changed = [column for column in PORTABLE_LOCATOR_COLUMNS if column in frame.columns]
    for column in changed:
        frame[column] = frame[column].map(portable_locator)
    return frame, changed


def export_snapshot(destination):
    """Copy the minimal redistributable inputs needed to rebuild the public release.

    All inputs are public metadata (GEO series matrices, PMC/publisher supplements,
    the authors' public metadata repository and the derived canonical tables);
    no expression matrix or controlled-access material is included.
    """
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    manifest = []
    for rel, origin in SNAPSHOT_INPUTS:
        source = PROJECT / rel
        if not source.exists():
            raise FileNotFoundError(source)
        target = destination / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        if source.suffix == '.parquet':
            frame, changed = portable_frame(source)
            frame.to_parquet(target, index=False)
        else:
            changed = []
            target.write_bytes(source.read_bytes())
        manifest.append(dict(path=rel, public_origin=origin, bytes=target.stat().st_size,
                             **({'locator_columns_converted': changed} if changed else {})))
    (destination / 'SNAPSHOT_MANIFEST.json').write_text(
        json.dumps(dict(purpose='minimal inputs to rebuild PAIR_LONGITUDINAL_MASTER via dataset/src/build_pair_dataset.py',
                        redistribution='public metadata only; no expression matrices, no controlled-access material',
                        files=manifest), indent=2) + '\n', encoding='utf-8')
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--release', default='v0.1.2')
    parser.add_argument('--export-snapshot', metavar='DIR',
                        help='copy the minimal redistributable inputs needed to rebuild the release')
    args = parser.parse_args()
    if args.export_snapshot:
        manifest = export_snapshot(args.export_snapshot)
        print(json.dumps(dict(snapshot_files=len(manifest),
                              snapshot_bytes=sum(f['bytes'] for f in manifest))))
        return 0
    archive = load_archive()
    tables, samples, endpoints, transitions, exposures, resolved_counts = resolved_frames(archive)
    patients = patient_index(archive, tables, samples, endpoints, transitions)
    master = build_master(archive, patients, samples, endpoints, transitions, exposures, tables['treatments'])
    dupes = duplicate_groups(patients)
    crosswalk = gse20181_gse5462_crosswalk(samples)
    prov = provenance_table(master, endpoints)
    summary = summarize(master, patients, samples, crosswalk)
    out = DATASET_ROOT / 'releases' / args.release
    out.mkdir(parents=True, exist_ok=True)
    # Written last and converted only here: everything above needs the real local path to
    # decide presence, and only the published tables must be portable.
    samples = samples.assign(expression_file=samples.expression_file.map(portable_locator))
    endpoints = endpoints.assign(source_file=endpoints.source_file.map(portable_locator))
    master.to_parquet(out / 'PAIR_LONGITUDINAL_MASTER.parquet', index=False)
    master.to_csv(out / 'PAIR_LONGITUDINAL_MASTER.csv.gz', index=False,
                  compression={'method': 'gzip', 'mtime': 0})
    patients.to_parquet(out / 'patients.parquet', index=False)
    samples.to_parquet(out / 'samples.parquet', index=False)
    endpoints.to_parquet(out / 'endpoints.parquet', index=False)
    prov.to_parquet(out / 'provenance.parquet', index=False)
    dupes.to_csv(out / 'duplicate_groups.csv', index=False)
    crosswalk.to_csv(out / 'gse20181_gse5462_crosswalk.csv', index=False)
    ledger_src = CORPUS / 'LABEL_RECOVERY_LEDGER.csv'
    if ledger_src.exists():
        (out / 'LABEL_RECOVERY_LEDGER.csv').write_bytes(ledger_src.read_bytes())
    (out / 'DATASET_SUMMARY.json').write_text(json.dumps(summary, indent=2, default=str) + '\n', encoding='utf-8')
    print(json.dumps({k: v for k, v in summary.items() if k != 'by_cohort'}, indent=2, default=str))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
