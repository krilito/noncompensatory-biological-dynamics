"""Targeted v0.1.2 invariants only — not a general test suite."""
import re
import sys
from pathlib import Path

import pandas as pd
import pytest

DATASET_ROOT = Path(__file__).resolve().parents[1]
RELEASE = DATASET_ROOT / 'releases' / 'v0.1.2'
sys.path.insert(0, str(DATASET_ROOT / 'src'))

import build_pair_dataset as bp  # noqa: E402

# The lookbehinds matter: without them the scheme separator of every https URL reads as a
# drive letter, and a URL is a portable locator, which is exactly what these columns should
# hold.  Every backslash and personal-directory root is assembled from fragments, because a
# literal copy of a scanner marker in this file would itself be flagged by
# scripts/scan_public_repository.py.
BACKSLASH = chr(92)
POSIX_ROOTS = [prefix + chr(47) for prefix in ('/ho' + 'me', '/U' + 'sers', '/mn' + 't')]
MACHINE_PATH = re.compile('|'.join([
    r'(?<![A-Za-z0-9._-])[A-Za-z][:]' + '[' + BACKSLASH * 2 + '/]',        # drive-qualified
    r'(?<![A-Za-z0-9._-])' + BACKSLASH * 4 + r'[A-Za-z0-9.$_-]',            # UNC share
    r'(?<![A-Za-z0-9._-])(?:' + '|'.join(POSIX_ROOTS) + ')']))             # personal roots


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


def _text_cells(path: Path):
    """Yield (column, value) for every text cell of a published tabular artifact."""
    if path.suffix == '.parquet':
        frame = pd.read_parquet(path)
    elif path.name.endswith('.csv.gz'):
        frame = pd.read_csv(path, dtype=str, low_memory=False)
    else:
        return
    for column in frame.columns:
        series = frame[column].dropna()
        if not (pd.api.types.is_object_dtype(series) or pd.api.types.is_string_dtype(series)):
            continue
        for value in series.astype(str):
            yield column, value


def test_no_published_dataset_artifact_carries_a_machine_absolute_path():
    """A release must locate a source by project-relative path or accession, never by drive.

    Machine location is not record identity.  The two columns that used to hold absolute
    paths (``samples.expression_file`` and ``endpoints.source_file``) are now written in
    project-relative POSIX form by ``portable_locator``, and an out-of-tree absolute path is
    labelled SOURCE_EXTERNAL instead of being invented.
    """
    offenders = []
    for path in sorted(DATASET_ROOT.rglob('*')):
        if not path.is_file() or '__pycache__' in path.parts:
            continue
        relative = path.relative_to(DATASET_ROOT).as_posix()
        if path.suffix == '.parquet' or path.name.endswith('.csv.gz'):
            found = [(relative, column, value[:90])
                     for column, value in _text_cells(path) if MACHINE_PATH.search(value)]
        elif path.suffix in {'.csv', '.json', '.md'}:
            text = path.read_text(encoding='utf-8', errors='replace')
            found = [(relative, '', match.group(0)) for match in
                     (MACHINE_PATH.search(text),) if match]
        else:
            continue
        offenders.extend(found)
    assert not offenders, offenders[:5]


def test_a_sentinel_never_becomes_an_invented_locator(tmp_path, monkeypatch):
    """Resolving a "no file" marker against the working directory fabricates a path.

    Before v0.3.2 the 214 MASTER rows without expression carried
    ``longitudinal-data/dataset/NOT_FOUND`` — a locator assembled from wherever the build was
    started, which is neither portable nor true.
    """
    project = tmp_path / 'project'
    monkeypatch.setattr(bp, 'PROJECT', project)
    monkeypatch.chdir(tmp_path)

    assert bp.expression_reference('NOT_FOUND') == 'NOT_AVAILABLE'
    assert bp.expression_reference('') == 'NOT_AVAILABLE'
    assert bp.expression_reference(None) == 'NOT_AVAILABLE'
    # an already-relative locator is project-relative; it is normalised, never re-resolved
    assert bp.portable_locator('longitudinal-data/data/raw/X/x.txt.gz') \
        == 'longitudinal-data/data/raw/X/x.txt.gz'
    inside = project / 'data' / 'raw' / 'x.gz'
    assert bp.portable_locator(str(inside)) == 'data/raw/x.gz'
    # an absolute path outside the project tree is labelled, never published as a machine path
    outside = tmp_path / 'elsewhere' / 'x.gz'
    assert bp.portable_locator(str(outside)) == 'SOURCE_EXTERNAL'
    assert bp.portable_locator('NOT_FOUND') == 'NOT_FOUND'
