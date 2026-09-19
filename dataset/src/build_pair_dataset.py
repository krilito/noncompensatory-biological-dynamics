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
# Same study, shared GSM sample accessions; resource-level duplicate, patients not merged.
SAME_STUDY_RESOURCE_GROUPS = {'DUP_GSE20181_GSE5462': ['GSE20181', 'GSE5462']}


def load_archive():
    spec = importlib.util.spec_from_file_location('patient_archive', CORPUS / 'src/patient_archive.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


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
        patients.loc[patients.patient_uid.eq(uid), 'duplicate_resource_group'] = uid
    return patients


def build_master(patients, samples, endpoints, transitions, exposures, treatments):
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

    def relpath(value):
        try:
            return str(Path(value).resolve().relative_to(PROJECT)).replace('\\', '/')
        except (ValueError, OSError):
            return str(value)

    patient_by_key = patients.set_index('patient_key')
    rows = []
    for tr in transitions.itertuples(index=False):
        key = tr.patient_key if isinstance(tr.patient_key, str) else ''
        patient = patient_by_key.loc[key] if key in patient_by_key.index else None
        s0, s1 = sample_index.loc[tr.sample_t0], sample_index.loc[tr.sample_t1]
        native_regimen = exposure_by_sample.get(tr.sample_t1, '')
        drug, drug_class, family = nz.classify_treatment(native_regimen, classes_by_sample.get(tr.sample_t1, ''))
        if patient is not None:
            status = patient.endpoint_status
            families = [f for f in str(patient.endpoint_families).split('|') if f]
            values = [v for v in str(patient.endpoint_native_values).split('|') if v]
            harmonized = str(patient.endpoint_harmonized)
            timing = str(patient.endpoint_assessment_time)
            verification = patient.endpoint_verification
            uid, identity, dupe = patient.patient_uid, patient.identity_status, patient.duplicate_resource_group
            cancer, cancer_native, subtype = patient.cancer_type, patient.cancer_type_native, patient.oncotree
        else:
            status, families, values, harmonized = 'MAPPING_NOT_FOUND', [], '', ''
            timing, verification = '', 'NOT_LINKED'
            uid, identity, dupe = key or 'UNMAPPED', 'MAPPING_NOT_FOUND', ''
            cancer = cancer_native = 'NOT_ESTABLISHED'
            subtype = 'NOT_MEASURED'
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
            phase_t0=nz.normalize_phase(s0.timepoint_normalized),
            phase_t1=nz.normalize_phase(s1.timepoint_normalized),
            delta_days='NOT_MEASURED',
            treatment_native=native_regimen or 'NOT_MEASURED',
            drug=drug, drug_class=drug_class, treatment_family=family,
            endpoint_family='|'.join(families),
            endpoint_system=nz.endpoint_systems(families),
            endpoint_native='|'.join(values),
            endpoint_harmonized=harmonized,
            endpoint_assessment_time=timing,
            endpoint_status=status,
            endpoint_verification_status=verification,
            expression_modality=s1.modality or 'NOT_MEASURED',
            platform=s1.platform or 'NOT_MEASURED',
            expression_t0_ref=relpath(s0.expression_file) if s0.expression_file else 'NOT_AVAILABLE',
            expression_t1_ref=relpath(s1.expression_file) if s1.expression_file else 'NOT_AVAILABLE',
            gene_space=nz.gene_space_for(s1.modality),
            source_accession=s1.dataset_id or tr.cohort_code,
            source_publication=nz.SOURCE_PUBLICATION.get(tr.cohort_code, ''),
            source_table='longitudinal-data/corpus/data/transitions.parquet',
            provenance_id=tr.transition_id,
            paired_valid=bool(tr.pair_valid),
            clinical_endpoint_available=status.startswith('OBSERVED_'),
            strict_auo_eligible=bool(tr.pair_valid) and status == 'OBSERVED_STRICT_RECIST',
        ))
    return pd.DataFrame(rows)


def duplicate_groups(patients):
    rows = []
    for gid, group in patients.groupby('duplicate_resource_group', sort=True):
        if not gid:
            continue
        members = sorted(group.patient_key)
        same_patient = group.patient_uid.nunique() == 1
        rows.append(dict(
            duplicate_group_id=gid,
            group_type='SAME_PATIENT' if same_patient else 'SAME_STUDY_RESOURCE',
            n_members=len(members), members='|'.join(members),
            evidence='GLOBAL_PATIENT_IDENTITY_CONFIRMED_SAME_PATIENT' if same_patient
            else 'SHARED_GSM_SAMPLE_ACCESSIONS_SAME_STUDY',
            patients_merged='YES' if same_patient else 'NO'))
    return pd.DataFrame(rows)


def provenance_table(master, endpoints):
    endpoint_sources = {}
    for key, group in endpoints[endpoints.patient_key.ne('')].groupby('patient_key'):
        endpoint_sources[key] = '|'.join(sorted(
            {Path(p).name for v in group.source_file for p in str(v).split('|') if p}))
    rows = []
    for row in master.itertuples(index=False):
        key = f'{row.cohort_code}:{row.native_patient_id}'
        families = set(row.endpoint_family.split('|'))
        rows.append(dict(
            pair_uid=row.pair_uid, patient_key=key, patient_uid=row.patient_uid,
            transition_source='corpus/data/transitions.parquet',
            endpoint_source_files=endpoint_sources.get(key, ''),
            endpoint_source='corpus/clinical/clinical_endpoints.parquet + recovery sources (see LABEL_RECOVERY_LEDGER.csv)',
            recovery_labels_used='YES' if families & RECOVERED_FAMILIES else 'NO'))
    return pd.DataFrame(rows)


def summarize(master, patients, samples):
    endpoint_verified = patients.endpoint_status.str.startswith('OBSERVED_')
    per_cohort = []
    for cohort, group in patients.groupby('cohort_code', sort=True):
        master_rows = master[master.cohort_code.eq(cohort)]
        per_cohort.append(dict(
            cohort_code=cohort, resource_patient_entries=len(group),
            unique_patients=group.patient_uid.nunique(),
            master_rows=len(master_rows),
            paired_patients=int(group.has_pair.sum()),
            endpoint_verified_patients=int(endpoint_verified[group.index].sum()),
            paired_endpoint_verified=int((group.has_pair & endpoint_verified[group.index]).sum()),
            strict_auo_eligible_rows=int(master_rows.strict_auo_eligible.sum()),
            cancer_type=group.cancer_type.iloc[0]))
    dupes = patients[patients.duplicate_resource_group.ne('')]
    return dict(
        pair_dataset_version='v0.1',
        generated_from='longitudinal-data/corpus canonical parquet tables + patient_archive recovery pipeline',
        cohort_count=len(per_cohort),
        resource_patient_entries=len(patients),
        unique_patients_after_confirmed_dedup=int(patients.patient_uid.nunique()),
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
        strict_auo_eligible_patients=int(master.groupby('patient_uid').strict_auo_eligible.any().sum()),
        frozen_auo_paper_subset=27,
        duplicate_resource_groups=int(dupes.duplicate_resource_group.nunique()),
        by_cohort=per_cohort,
        note='unique_patients reflects only confirmed cross-resource identities (GIDE/MORRISON_gide, '
             'GSE91061/MORRISON_038); GSE20181/GSE5462 is a same-study resource duplicate whose '
             'patients are not merged')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--release', default='v0.1')
    args = parser.parse_args()
    archive = load_archive()
    tables, samples, endpoints, transitions, exposures, resolved_counts = resolved_frames(archive)
    patients = patient_index(archive, tables, samples, endpoints, transitions)
    master = build_master(patients, samples, endpoints, transitions, exposures, tables['treatments'])
    dupes = duplicate_groups(patients)
    prov = provenance_table(master, endpoints)
    summary = summarize(master, patients, samples)
    out = DATASET_ROOT / 'releases' / args.release
    out.mkdir(parents=True, exist_ok=True)
    master.to_parquet(out / 'PAIR_LONGITUDINAL_MASTER.parquet', index=False)
    master.to_csv(out / 'PAIR_LONGITUDINAL_MASTER.csv.gz', index=False,
                  compression={'method': 'gzip', 'mtime': 0})
    patients.to_parquet(out / 'patients.parquet', index=False)
    samples.to_parquet(out / 'samples.parquet', index=False)
    endpoints.to_parquet(out / 'endpoints.parquet', index=False)
    prov.to_parquet(out / 'provenance.parquet', index=False)
    dupes.to_csv(out / 'duplicate_groups.csv', index=False)
    ledger_src = CORPUS / 'LABEL_RECOVERY_LEDGER.csv'
    if ledger_src.exists():
        (out / 'LABEL_RECOVERY_LEDGER.csv').write_bytes(ledger_src.read_bytes())
    (out / 'DATASET_SUMMARY.json').write_text(json.dumps(summary, indent=2, default=str) + '\n', encoding='utf-8')
    print(json.dumps({k: v for k, v in summary.items() if k != 'by_cohort'}, indent=2, default=str))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
