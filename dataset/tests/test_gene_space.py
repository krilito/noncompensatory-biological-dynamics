"""v0.2 canonical gene-space invariants — artifact level, not a general test suite."""
import json
from pathlib import Path

import pandas as pd
import pytest

GENE_SPACE = Path(__file__).resolve().parents[1] / 'gene_space'
VOCABULARY = {'EXACT', 'ALIAS', 'PREVIOUS_SYMBOL', 'PLATFORM_ANNOTATION', 'AMBIGUOUS', 'UNMAPPED'}
MAPPED = VOCABULARY - {'AMBIGUOUS', 'UNMAPPED'}
FEATURE_COLUMNS = ['cohort_code', 'source_expression', 'platform', 'original_feature_id',
                   'original_feature_type', 'original_feature_annotation',
                   'feature_identifier_class', 'gene_addressable', 'canonical_gene_id',
                   'hgnc_symbol', 'ensembl_gene_id', 'entrez_gene_id', 'mapping_method',
                   'mapping_status', 'mapping_ambiguity', 'reference_version']


def report():
    path = GENE_SPACE / 'GENE_SPACE_REPORT.json'
    if not path.exists():
        pytest.skip('gene space not built')
    return json.loads(path.read_text(encoding='utf-8'))


def metrics():
    return pd.read_csv(GENE_SPACE / 'SOURCE_MAPPING_METRICS.csv')


def feature_map(source):
    return pd.read_csv(GENE_SPACE / 'feature_map' / f'{source}.feature_map.csv.gz', dtype=str)


def mapped_sources():
    return metrics().query("status == 'MAPPED'").source_expression


def test_every_mapped_source_has_a_feature_map_with_the_canonical_schema():
    for source in mapped_sources():
        path = GENE_SPACE / 'feature_map' / f'{source}.feature_map.csv.gz'
        assert path.exists(), source
        assert list(pd.read_csv(path, nrows=1).columns) == FEATURE_COLUMNS, source


def test_status_vocabulary_is_closed_and_ambiguity_is_never_resolved_by_choice():
    for source in mapped_sources():
        frame = feature_map(source)
        assert set(frame.mapping_status) <= VOCABULARY, source
        ambiguous = frame[frame.mapping_status.eq('AMBIGUOUS')]
        assert ambiguous.hgnc_symbol.isna().all(), source
        assert ambiguous.canonical_gene_id.isna().all(), source
        assert ambiguous.mapping_ambiguity.notna().all(), source
        assert (ambiguous.mapping_ambiguity.str.count(r'\|') >= 1).all(), source
        mapped = frame[frame.mapping_status.isin(MAPPED)]
        assert mapped.hgnc_symbol.notna().all(), source
        assert mapped.canonical_gene_id.str.startswith('HGNC:').all(), source
        assert not mapped.hgnc_symbol.str.contains(r'\|').any(), source


def test_native_rows_are_preserved_even_when_the_class_is_not_gene_addressable():
    for source in mapped_sources():
        frame = feature_map(source)
        unit = metrics().query('source_expression == @source').iloc[0]
        assert len(frame) == unit.raw_feature_count, source
        assert frame.original_feature_id.notna().all(), source
        assert set(frame.source_expression) == {source}
        assert frame.reference_version.notna().all(), source
        excluded = frame[frame.gene_addressable.eq('False')]
        assert len(excluded) == unit.non_addressable_feature_count, source
        # An excluded row is still published with its native identifier and its own verdict.
        assert excluded.original_feature_id.ne('').all(), source
        assert excluded.feature_identifier_class.ne('GENE_IDENTIFIER').all(), source
        assert set(excluded.feature_identifier_class) <= set(json.loads(
            unit.non_addressable_classes)), source


def test_both_denominators_are_reported_and_never_conflate():
    for _, unit in metrics().query("status == 'MAPPED'").iterrows():
        counts = feature_map(unit.source_expression).mapping_status.value_counts()
        assert unit.n_unmapped == int(counts.get('UNMAPPED', 0)), unit.source_expression
        assert unit.raw_mapped_feature_count == unit.raw_feature_count - unit.n_ambiguous \
            - unit.n_unmapped, unit.source_expression
        assert abs(unit.raw_mapping_fraction
                   - unit.raw_mapped_feature_count / unit.raw_feature_count) < 1e-4
        assert unit.gene_addressable_feature_count <= unit.raw_feature_count
        assert unit.gene_addressable_mapped_count <= unit.gene_addressable_feature_count
        # The gene-addressable numerator can only contain gene-addressable features.
        assert unit.gene_addressable_mapped_count <= unit.raw_mapped_feature_count, \
            unit.source_expression
        assert unit.exact_mapping_fraction <= unit.raw_mapping_fraction + 1e-9
        assert unit.raw_feature_count - unit.gene_addressable_feature_count \
            == unit.non_addressable_feature_count, unit.source_expression


def test_core_entry_is_decided_by_criteria_not_by_source_names():
    frame = metrics().query("status == 'MAPPED'")
    expected = frame.raw_coverage_passes_entry_rule | (
        frame.gene_addressable_coverage_passes_entry_rule & frame.addressable_admission_evidenced)
    assert set(frame.loc[expected, 'source_expression']) == \
        set(frame.loc[frame.enters_core_gene_space, 'source_expression'])
    # Admission on the addressable denominator always carries namespace evidence.
    admitted = frame[frame.enters_core_gene_space & ~frame.raw_coverage_passes_entry_rule]
    assert admitted.addressable_admission_evidenced.all()
    # Missing platform annotation is not an excuse: GSE3578 stays below the rule.
    gse3578 = frame.query("source_expression == 'GSE3578_array'").iloc[0]
    assert not gse3578.enters_core_gene_space
    assert not gse3578.raw_coverage_passes_entry_rule
    evidence = json.loads(gse3578.non_addressable_classes)
    assert not any(ev['outside_hgnc_gene_namespace']
                   for name, ev in evidence.items() if name == 'unannotated_array_feature')


def test_low_raw_coverage_admissions_are_evidenced():
    report_frame = report()
    admitted = report_frame['adjudication']['admitted_by_gene_addressable_rule']
    for source in admitted:
        unit = metrics().query('source_expression == @source').iloc[0]
        classes = json.loads(unit.non_addressable_classes)
        outside = {name: ev for name, ev in classes.items() if ev['outside_hgnc_gene_namespace']}
        assert outside, source
        # Evidence: the excluded classes are near-absent from the HGNC namespace, and their
        # feature count is exactly what shrank the denominator.
        assert all(ev['namespace_presence_rate'] <= 0.02 for ev in outside.values()), source
        assert sum(ev['features'] for ev in outside.values()) \
            == unit.non_addressable_feature_count, source
        assert unit.raw_mapping_fraction < 0.5 <= unit.gene_addressable_mapping_fraction, source


def test_core_gene_space_is_a_subset_of_the_union_and_agrees_with_the_report():
    genes = pd.read_csv(GENE_SPACE / 'GENE_SPACE_GENES.csv.gz')
    headline = report()['headline']
    assert len(genes) == headline['genes_in_union']
    assert int(genes.in_core_gene_space.sum()) == headline['genes_in_core_all_eligible_sources']
    n_eligible = headline['sources_entering_core']
    assert genes.loc[genes.in_core_gene_space, 'n_eligible_sources'].eq(n_eligible).all()
    assert genes.loc[~genes.in_core_gene_space, 'n_eligible_sources'].lt(n_eligible).all()
    assert genes.hgnc_symbol.is_unique
    # The core is not the naive all-platform intersection.
    assert headline['genes_in_core_all_eligible_sources'] > headline['genes_in_intersection_all_sources']


def test_pair_gene_space_counts_are_reproducible_from_the_artifact():
    pairs = pd.read_csv(GENE_SPACE / 'PAIR_GENE_SPACE.csv.gz')
    headline = report()['headline']
    with_expression = pairs[~pairs.gene_space_status.eq('EXPRESSION_FILE_NOT_LOCAL')]
    assert len(with_expression) == headline['pairs_with_t0_t1_expression']
    assert int(pairs.in_canonical_gene_space.sum()) == headline['pairs_with_canonical_gene_space']
    core = set(metrics().query("enters_core_gene_space == True").source_expression)
    assert set(with_expression.loc[pairs.in_canonical_gene_space, 'source_expression']) <= core
    excluded = with_expression[~with_expression.in_canonical_gene_space]
    assert not excluded.source_expression.isin(core).any()
    # Every excluded interval is named by a rule outcome, never silently dropped.
    assert excluded.gene_space_status.str.startswith(('BELOW_CORE_ENTRY_RULE', 'MIXED',
                                                       'UNDECLARED')).all()


def test_no_machine_bound_path_reaches_a_gene_space_text_artifact():
    # No artifact may carry a Windows path separator or a drive-qualified path.
    for path in list(GENE_SPACE.glob('*.csv')) + list(GENE_SPACE.glob('*.json')):
        assert chr(92) not in path.read_text(encoding='utf-8'), path.name
    # Every path written into an artifact is repository-relative and POSIX-separated.
    for source in metrics().query("status == 'MAPPED'").source_expression:
        file = metrics().query('source_expression == @source').iloc[0].expression_file
        assert chr(92) not in file and ':' not in file[:3], source


def test_gene_space_artifacts_pass_the_public_repository_patterns():
    """In a published copy, run the public tree's own fail-closed patterns over the artifacts."""
    scripts = Path(__file__).resolve().parents[2] / 'scripts'
    if not (scripts / 'scan_public_repository.py').exists():
        pytest.skip('not running from a published copy')
    import sys
    sys.path.insert(0, str(scripts))
    from scan_public_repository import ABSOLUTE_PATH_RE, PRIVATE_IDENTITY_RE, SECRET_RE
    paths = (list(GENE_SPACE.glob('*.csv')) + list(GENE_SPACE.glob('*.json'))
             + list(GENE_SPACE.glob('*.md')) + [Path(__file__)])
    for path in paths:
        text = path.read_text(encoding='utf-8')
        for pattern in (ABSOLUTE_PATH_RE, PRIVATE_IDENTITY_RE, SECRET_RE):
            assert not pattern.search(text), (path.name, pattern.pattern)
