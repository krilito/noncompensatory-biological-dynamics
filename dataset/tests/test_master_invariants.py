"""Targeted v0.1.2 invariants only — not a general test suite."""
from pathlib import Path

import pandas as pd
import pytest

RELEASE = Path(__file__).resolve().parents[1] / 'releases' / 'v0.1.2'


def master():
    path = RELEASE / 'PAIR_LONGITUDINAL_MASTER.parquet'
    if not path.exists():
        pytest.skip('v0.1.2 release not built')
    return pd.read_parquet(path)


def test_endpoints_bind_to_intervals_not_copied():
    m = master()
    assert set(m.columns) >= {'patient_endpoint_status', 'interval_endpoint_status'}
    multi = m[m.duplicated('patient_uid', keep=False)]
    unresolved = multi[multi.interval_endpoint_status.eq('NOT_INTERVAL_RESOLVED')]
    # Unresolved intervals must not carry interval endpoint values.
    assert (unresolved.interval_endpoint_native.fillna('').eq('')).all()
    # Patient-level outcome is never dropped by interval binding.
    assert (multi.patient_endpoint_status.ne('LABEL_NOT_FOUND') |
            multi.patient_endpoint_native.fillna('').ne('')).any()


def test_frozen_auo_is_27_gse91061_only():
    m = master()
    frozen = m[m.frozen_auo_eligible]
    assert frozen.patient_uid.nunique() == 27
    assert set(frozen.cohort_code) == {'GSE91061'}


def test_cross_dataset_eligibility_is_a_separate_field():
    m = master()
    assert 'strict_prcr_vs_pd_eligible' in m.columns
    assert 'strict_auo_eligible' not in m.columns
    strict = m[m.strict_prcr_vs_pd_eligible]
    assert strict.patient_uid.nunique() > 27  # broader than the frozen subset
    assert set(strict.cohort_code) > {'GSE91061'}


def test_gse20181_gse5462_patients_merged_at_uid_level():
    m = master()
    import json
    summary = json.loads((RELEASE / 'DATASET_SUMMARY.json').read_text(encoding='utf-8'))
    assert summary['gsm_confirmed_merged_patient_pairs'] == 58
    group = summary['same_study_duplicate_resources'][0]
    assert group['duplicate_group_id'] == 'DUP_GSE20181_GSE5462'
    assert group['patients_merged'] == 'YES'
    # Every shared-GSM patient pair carries one patient_uid across both resources.
    crosswalk = pd.read_csv(RELEASE / 'gse20181_gse5462_crosswalk.csv', dtype=str)
    shared = crosswalk[crosswalk.evidence.eq('SHARED_GSM_SAMPLE_ACCESSION')]
    assert len(shared) == 116
    for row in shared.drop_duplicates(['patient_gse20181', 'patient_gse5462']).itertuples(index=False):
        left = m[m.cohort_code.eq('GSE20181') & m.native_patient_id.eq(str(row.patient_gse20181))]
        right = m[m.cohort_code.eq('GSE5462') & m.native_patient_id.eq(str(row.patient_gse5462))]
        assert left.patient_uid.nunique() == 1 and right.patient_uid.nunique() == 1
        assert left.patient_uid.iloc[0] == right.patient_uid.iloc[0]
