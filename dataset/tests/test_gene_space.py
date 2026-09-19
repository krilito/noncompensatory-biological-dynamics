"""v0.2 canonical gene-space invariants — artifact level, not a general test suite."""
import json
from pathlib import Path

import pandas as pd
import pytest

GENE_SPACE = Path(__file__).resolve().parents[1] / 'gene_space'
VOCABULARY = {'EXACT', 'ALIAS', 'PREVIOUS_SYMBOL', 'PLATFORM_ANNOTATION', 'AMBIGUOUS', 'UNMAPPED'}
MAPPED = VOCABULARY - {'AMBIGUOUS', 'UNMAPPED'}


def report():
    path = GENE_SPACE / 'GENE_SPACE_REPORT.json'
    if not path.exists():
        pytest.skip('gene space not built')
    return json.loads(path.read_text(encoding='utf-8'))


def metrics():
    return pd.read_csv(GENE_SPACE / 'SOURCE_MAPPING_METRICS.csv')


def feature_map(source):
    return pd.read_csv(GENE_SPACE / 'feature_map' / f'{source}.feature_map.csv.gz', dtype=str)


def test_every_mapped_source_has_a_feature_map_with_the_canonical_schema():
    columns = ['cohort_code', 'source_expression', 'platform', 'original_feature_id',
               'original_feature_type', 'original_feature_annotation', 'canonical_gene_id',
               'hgnc_symbol', 'ensembl_gene_id', 'entrez_gene_id', 'mapping_method',
               'mapping_status', 'mapping_ambiguity', 'reference_version']
    for source in metrics().query("status == 'MAPPED'").source_expression:
        path = GENE_SPACE / 'feature_map' / f'{source}.feature_map.csv.gz'
        assert path.exists(), source
        assert list(pd.read_csv(path, nrows=1).columns) == columns, source


def test_status_vocabulary_is_closed_and_ambiguity_is_never_resolved_by_choice():
    for source in metrics().query("status == 'MAPPED'").source_expression:
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


def test_native_representation_is_preserved_for_every_feature():
    for source in metrics().query("status == 'MAPPED'").source_expression:
        frame = feature_map(source)
        unit = metrics().query('source_expression == @source').iloc[0]
        assert len(frame) == unit.source_features, source
        assert frame.original_feature_id.notna().all(), source
        assert set(frame.source_expression) == {source}
        assert frame.reference_version.notna().all(), source


def test_per_source_counts_match_the_report_fractions():
    for _, unit in metrics().query("status == 'MAPPED'").iterrows():
        counts = feature_map(unit.source_expression).mapping_status.value_counts()
        assert unit.n_ambiguous == int(counts.get('AMBIGUOUS', 0)), unit.source_expression
        assert unit.n_unmapped == int(counts.get('UNMAPPED', 0)), unit.source_expression
        assert unit.mapped_features == unit.source_features - unit.n_ambiguous - unit.n_unmapped
        assert abs(unit.mapped_fraction - unit.mapped_features / unit.source_features) < 1e-4
        assert unit.exact_mapping_fraction <= unit.mapped_fraction + 1e-9


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
    # A pair never enters the canonical space through a source that failed the entry rule.
    excluded = with_expression[~with_expression.in_canonical_gene_space]
    assert not excluded.source_expression.isin(core).any()


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
