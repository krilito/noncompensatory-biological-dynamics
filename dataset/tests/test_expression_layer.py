"""Level C1 checks: the decisions that would silently corrupt the data.

These tests are deliberately small.  They pin the *order* of aggregation versus
library-size normalization, the native origin of the CPM denominator, the collapse rule
for duplicate mapped features, the no-double-log guarantee, the fail-closed behaviour for
undetermined scales, that a sample naming several matrix columns is never collapsed, and
the internal consistency of the published pair-coverage artifact.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))

import build_expression_layer as builder
from expression_transforms import (
    ExpressionScale, SourceExpressionContract, build_canonical_quantitative_matrix)

LAYER = builder.LAYER

FEATURE_COLUMNS = ['original_feature_id', 'hgnc_symbol', 'mapping_status']


def contract(scale: ExpressionScale, source: str = 'TEST') -> SourceExpressionContract:
    return SourceExpressionContract(
        source_expression=source, modality='test', scale=scale,
        gene_aggregation=('SUM_BEFORE_TRANSFORM' if scale is ExpressionScale.RAW_COUNTS
                          else 'MEDIAN_AFTER_TRANSFORM'))


def frame(values: dict[str, list[float]], index: list[str]) -> pd.DataFrame:
    return pd.DataFrame(values, index=index)


def feature_map(rows: list[tuple[str, str]]) -> pd.DataFrame:
    return pd.DataFrame(rows, columns=FEATURE_COLUMNS)


def test_counts_are_summed_to_gene_before_cpm():
    """Summing after CPM would change the additive structure the counts carry."""
    native = frame({'sample_a': [10.0, 30.0, 60.0], 'sample_b': [20.0, 40.0, 80.0]},
                   ['F1', 'F2', 'F3'])
    mapping = feature_map([('F1', 'GENE_X', 'EXACT'), ('F2', 'GENE_X', 'EXACT'),
                           ('F3', 'GENE_Y', 'EXACT')])
    output, report = build_canonical_quantitative_matrix(
        native, mapping, contract(ExpressionScale.RAW_COUNTS))

    # sample_a library size is 100 and GENE_X is 10 + 30 counts, not (0.1 + 0.3) CPM.
    assert report['transform'] == 'SUM_TO_GENE_THEN_LOG2_NATIVE_CPM_PLUS_1'
    assert output.loc['GENE_X', 'sample_a'] == pytest.approx(np.log2(40.0 / 100.0 * 1e6 + 1))
    assert output.loc['GENE_X', 'sample_b'] == pytest.approx(np.log2(60.0 / 140.0 * 1e6 + 1))
    assert output.loc['GENE_Y', 'sample_a'] == pytest.approx(np.log2(600_000.0 + 1))


def test_unmapped_native_features_still_consume_library_depth():
    """The v0.3.1 fix: the CPM denominator is the native matrix, not what survived mapping.

    A feature with no canonical gene is absent from the output but its reads were really
    sequenced.  Dividing by the mapped mass only would rescale every retained gene by the
    inverse of the mapping fraction, silently renormalizing the source we were told not to
    renormalize.
    """
    native = frame({'sample_a': [10.0, 30.0, 5.0, 55.0]}, ['F1', 'F2', 'F3', 'NONHSAG1'])
    mapping = feature_map([('F1', 'GENE_X', 'EXACT'), ('F2', 'GENE_X', 'EXACT'),
                           ('F3', 'GENE_Y', 'EXACT'),
                           ('NONHSAG1', '', 'UNMAPPED')])
    output, report = build_canonical_quantitative_matrix(
        native, mapping, contract(ExpressionScale.RAW_COUNTS))

    assert 'NONHSAG1' not in output.index
    # native mass is 100, mapped mass is 45: GENE_X is 40/100, not 40/45.
    assert output.loc['GENE_X', 'sample_a'] == pytest.approx(np.log2(400_000.0 + 1))
    assert report['native_count_mass']['sample_a'] == 100.0
    assert report['mapped_count_mass']['sample_a'] == 45.0
    linear = np.power(2.0, output['sample_a'].to_numpy(dtype=float)) - 1.0
    assert linear.sum() == pytest.approx(1e6 * 45.0 / 100.0)
    assert linear.sum() != pytest.approx(1e6)


def test_array_probes_are_median_collapsed_and_never_summed():
    """Probe intensity is not additive: summing three probes would triple the gene."""
    native = frame({'sample_a': [5.0, 7.0, 9.0], 'sample_b': [1.0, 2.0, 3.0]},
                   ['PROBE1', 'PROBE2', 'PROBE3'])
    mapping = feature_map([(probe, 'TP53', 'PLATFORM_ANNOTATION')
                           for probe in ['PROBE1', 'PROBE2', 'PROBE3']])
    output, _ = build_canonical_quantitative_matrix(
        native, mapping, contract(ExpressionScale.ARRAY_NORMALIZED_LOG))

    assert output.loc['TP53', 'sample_a'] == pytest.approx(7.0)
    assert output.loc['TP53', 'sample_a'] != pytest.approx(21.0)
    assert output.loc['TP53', 'sample_b'] == pytest.approx(2.0)


def test_already_transformed_scales_are_not_logged_a_second_time():
    """The single most damaging possible bug: logarithming an author-supplied logCPM."""
    for scale in (ExpressionScale.LOG2_CPM, ExpressionScale.LOG2_TPM,
                  ExpressionScale.LOG2_FPKM, ExpressionScale.LOG_EXPRESSION,
                  ExpressionScale.ARRAY_NORMALIZED_LOG):
        native = frame({'sample_a': [-4.2, 120.0], 'sample_b': [0.0, 7.0]},
                       ['G1', 'G2'])
        mapping = feature_map([('G1', 'GENE_A', 'EXACT'), ('G2', 'GENE_B', 'EXACT')])
        output, report = build_canonical_quantitative_matrix(native, mapping, contract(scale))
        assert report['transform'] == 'IDENTITY_THEN_MEDIAN_TO_GENE'
        assert output.loc['GENE_A', 'sample_a'] == pytest.approx(-4.2)
        assert output.loc['GENE_B', 'sample_b'] == pytest.approx(7.0)


def test_unknown_scale_fails_closed_instead_of_being_guessed():
    native = frame({'sample_a': [1.0, 2.0]}, ['G1', 'G2'])
    mapping = feature_map([('G1', 'GENE_A', 'EXACT'), ('G2', 'GENE_B', 'EXACT')])
    with pytest.raises(ValueError, match='UNKNOWN'):
        build_canonical_quantitative_matrix(native, mapping, contract(ExpressionScale.UNKNOWN))


def test_declared_contract_is_never_silently_repaired():
    """A contradiction between declaration and data is a failure, not a cue to re-declare."""
    continuous = frame({'sample_a': [1.5, 2.5]}, ['G1', 'G2'])
    mapping = feature_map([('G1', 'GENE_A', 'EXACT'), ('G2', 'GENE_B', 'EXACT')])
    with pytest.raises(ValueError, match='integer-like'):
        build_canonical_quantitative_matrix(continuous, mapping,
                                           contract(ExpressionScale.RAW_COUNTS))
    below_floor = frame({'sample_a': [-1.4, 2.5]}, ['G1', 'G2'])
    with pytest.raises(ValueError, match='floor'):
        build_canonical_quantitative_matrix(
            below_floor, mapping, contract(ExpressionScale.LINEAR_ABUNDANCE_COMBAT_ADJUSTED))


def test_published_pair_coverage_agrees_with_the_binding_manifest():
    coverage = pd.read_csv(LAYER / 'PAIR_EXPRESSION_COVERAGE.csv.gz')
    binding = pd.read_csv(LAYER / 'manifests/SAMPLE_COLUMN_BINDING.csv.gz')
    metrics = pd.read_csv(LAYER / 'SOURCE_EXPRESSION_METRICS.csv')
    report = json.loads((LAYER / 'EXPRESSION_LAYER_REPORT.json').read_text(encoding='utf-8'))

    bound = binding[binding.binding_status.eq('BOUND')]
    ready = set(metrics.loc[metrics.status.eq('QUANTITATIVE_READY'), 'source_expression'])
    usable = coverage.quantitative_status.eq('QUANTITATIVE_READY')

    assert int(usable.sum()) == report['headline']['pairs_quantitative_ready']
    assert coverage.loc[usable, 'source_expression'].isin(ready).all()
    # every ready pair has both endpoints bound, to two different matrix columns
    for row in coverage[usable].itertuples():
        assert row.matrix_column_t0 and row.matrix_column_t1
        assert row.matrix_column_t0 != row.matrix_column_t1
        assert {row.sample_t0, row.sample_t1}.issubset(set(bound.sample_id))
    # nothing is double counted and no pair is left with an unexplained status
    assert coverage.pair_uid.is_unique
    assert set(coverage.quantitative_status) <= {
        'QUANTITATIVE_READY', 'SEMANTICS_NOT_ESTABLISHED', 'SCALE_CONTRACT_FAILED',
        'SAMPLE_BINDING_FAILED', 'ROUTE_CHECK_FAILED', 'LIBRARY_DENOMINATOR_FAILED',
        'SOURCE_NOT_RESOLVED_IN_GENE_SPACE', 'ENDPOINT_SAMPLE_NOT_BOUND_TO_MATRIX_COLUMN'}


# HAZARDS is built from character codes: a literal copy of any marker in this file would
# itself be flagged by the public repository scanner.
HAZARDS = ['D' + chr(58) + chr(92), chr(47) + 'home' + chr(47),
           chr(47) + 'Users' + chr(47), chr(47) + 'mnt' + chr(47)]


def test_every_published_matrix_matches_its_declared_route_exactly():
    """ROUTE_IMPLEMENTATION_CHECK: was anything normalized or logarithmed twice?

    The builder recomputes each declared route from the native matrix and compares it
    with what was written, on genes served by exactly one mapped feature.  This proves
    the code followed the route; it does not prove the route is what the data deserves,
    which is the separate denominator test below.
    """
    metrics = pd.read_csv(LAYER / 'SOURCE_EXPRESSION_METRICS.csv')
    ready = metrics[metrics.status.eq('QUANTITATIVE_READY')]
    assert (ready.route_check_max_abs_diff == 0.0).all()
    assert (ready.route_check_genes > 500).all()
    assert set(ready.route_check_expected) == {
        'LOG2_OF_NATIVE_CPM_PLUS_1', 'LOG2_OF_X_PLUS_1', 'IDENTITY'}


def test_counts_matrices_keep_their_native_library_denominator():
    """LIBRARY_DENOMINATOR_CHECK: the canonical genes carry their true share, not 1e6.

    Reverting to a mapped-only denominator would make every counts column sum to exactly
    1e6 again, which is precisely the artifact v0.3 published as a strong verification.
    """
    metrics = pd.read_csv(LAYER / 'SOURCE_EXPRESSION_METRICS.csv')
    counts = metrics[metrics.input_scale.eq('RAW_COUNTS')]
    assert counts.status.eq('QUANTITATIVE_READY').all()
    assert counts.canonical_cpm_mass_identity_holds.notna().all()
    assert counts.canonical_cpm_mass_identity_holds.astype(bool).all()
    assert (counts.library_denominator_max_rel_dev < 1e-9).all()
    fractions = counts[['canonical_count_mass_fraction_min',
                        'canonical_count_mass_fraction_median',
                        'canonical_count_mass_fraction_max']].to_numpy(dtype=float)
    assert ((fractions > 0.0) & (fractions <= 1.0)).all()
    # the correction must be visible: some native features never entered the gene space
    assert (counts.canonical_count_mass_fraction_median < 0.999).any()
    assert (counts.linear_library_sum_max <= 1e6 + 1.0).all()
    # a source we were not entitled to renormalize must not look per-million.  The sums are
    # published for every ready source, so an empty column cannot pass this check.
    ready = metrics[metrics.status.eq('QUANTITATIVE_READY')]
    assert ready.linear_library_sum_min.notna().all()
    per_million = ready[ready.input_scale.ne('RAW_COUNTS')].linear_library_sum_min
    assert not per_million.map(lambda v: abs(v - 1e6) < 1.0).any()
    denominators = pd.read_csv(LAYER / 'manifests/RAW_COUNT_LIBRARY_DENOMINATORS.csv.gz')
    assert set(denominators.source_expression) == set(counts.source_expression)
    assert (denominators.mapped_count_mass <= denominators.native_count_mass).all()
    assert (denominators.native_count_mass > 0).all()


def test_a_sample_naming_several_matrix_columns_is_never_collapsed():
    """Multi-column rows are region, aliquot or unresolved cases, not libraries to add up."""
    binding = pd.read_csv(LAYER / 'manifests/SAMPLE_COLUMN_BINDING.csv.gz')
    adjudication = pd.read_csv(LAYER / 'manifests/SAMPLE_MULTI_COLUMN_ADJUDICATION.csv.gz')
    ambiguous = binding[binding.binding_status.eq('AMBIGUOUS_MULTIPLE_MATRIX_COLUMNS')]

    assert set(adjudication.sample_id) == set(ambiguous.sample_id)
    assert adjudication.sample_id.is_unique
    assert set(adjudication.replicate_type) <= set(builder.REPLICATE_TYPE_VOCABULARY)
    # nothing in C1 collapses: no summing, no averaging, and no bound multi-column row
    assert not adjudication.recommended_action.str.contains('SUM').any()
    assert not adjudication.replicate_type.eq('TECHNICAL_LIBRARY_REPLICATE').any()
    assert not ambiguous.sample_id.isin(binding[binding.binding_status.eq('BOUND')].sample_id).any()
    # a biological multi-region source must stay a region-level question, not a mean
    tracer = adjudication[adjudication.source_expression.eq('GSE139533_bulk')]
    assert set(tracer.replicate_type) == {'BIOLOGICAL_REGION_REPLICATE'}
    assert (tracer.n_candidate_columns > 1).all()


def test_author_corrected_source_is_flagged_for_the_model_layer():
    """C1 may keep an author-batch-corrected matrix, but must not hide what it is."""
    metrics = pd.read_csv(LAYER / 'SOURCE_EXPRESSION_METRICS.csv')
    combat = metrics[metrics.source_expression.eq('GSE319641_bulk')].iloc[0]
    assert bool(combat.author_batch_corrected)
    assert combat.preprocessing_scope == 'AUTHOR_FULL_SOURCE'
    assert combat.fold_isolation == 'NOT_ESTABLISHED'
    assert set(metrics.preprocessing_scope) <= {'OUR_C1_PER_SAMPLE', 'AUTHOR_FULL_SOURCE'}
    assert set(metrics.loc[metrics.preprocessing_scope.eq('OUR_C1_PER_SAMPLE'),
                           'fold_isolation']) == {'FOLD_INDEPENDENT'}


def test_series_matrix_lookup_is_deterministic_and_fails_closed(tmp_path, monkeypatch):
    """A documenting matrix is found by name, never guessed between two, never silently absent.

    Some submitters distribute per-platform matrices (GSE115821-GPL18573_...), so the plain
    accession name is not the only legal spelling; but returning {} used to make an absent
    matrix look like a source with no GEO labels, which silently binds zero samples.
    """
    project = tmp_path / 'project'

    def arrange(accession: str, names: list[str]) -> None:
        root = project / 'longitudinal-data' / 'data' / 'raw' / accession
        root.mkdir(parents=True, exist_ok=True)
        for name in names:
            (root / name).write_bytes(b'x')

    monkeypatch.setattr(builder.bg, 'PROJECT', project)
    monkeypatch.setattr(builder.bg, 'RAW', project / 'longitudinal-data' / 'data' / 'raw')

    arrange('GSE900001', ['GSE900001_series_matrix.txt.gz'])
    assert builder.resolve_series_matrix('GSE900001').name == 'GSE900001_series_matrix.txt.gz'

    arrange('GSE900002', ['GSE900002-GPL18573_series_matrix.txt.gz'])
    assert builder.resolve_series_matrix('GSE900002').name == 'GSE900002-GPL18573_series_matrix.txt.gz'

    arrange('GSE900003', ['GSE900003_series_matrix.txt.gz',
                          'GSE900003-GPL115821_series_matrix.txt.gz'])
    with pytest.raises(RuntimeError, match='AMBIGUOUS_DOCUMENTING_SERIES_MATRIX'):
        builder.resolve_series_matrix('GSE900003')

    arrange('GSE900004', [])
    with pytest.raises(FileNotFoundError, match='MISSING_DOCUMENTING_SERIES_MATRIX'):
        builder.resolve_series_matrix('GSE900004')


def test_a_source_that_does_not_read_geo_labels_never_probes_the_filesystem(monkeypatch):
    """The MGH sources bind by library label; an unrelated missing matrix must not stop them."""
    monkeypatch.setattr(builder.bg, 'RAW', Path('nowhere-at-all'))
    monkeypatch.setattr(builder.bg, 'PROJECT', Path('nowhere-at-all'))
    assert builder.geo_sample_labels('MGH_GSE115821_bulk') == {}
    assert builder.geo_sample_labels('MORRISON_bulk') == {}


def test_no_committed_artifact_carries_a_machine_bound_path():
    for path in sorted(LAYER.rglob('*')):
        if path.is_file() and path.suffix in {'.csv', '.json', '.md'}:
            text = path.read_text(encoding='utf-8')
            for hazard in HAZARDS:  # machine-bound markers, see below
                assert hazard not in text, (path, hazard)

def test_matrix_files_are_local_and_not_part_of_the_published_layer():
    published = {path.name for path in LAYER.iterdir()}
    assert published == {'SOURCE_EXPRESSION_CONTRACT.csv', 'SOURCE_EXPRESSION_METRICS.csv',
                         'EXPRESSION_LAYER_REPORT.json', 'PAIR_EXPRESSION_COVERAGE.csv.gz',
                         'manifests', 'README.md'}
    metrics = pd.read_csv(LAYER / 'SOURCE_EXPRESSION_METRICS.csv')
    assert set(builder.METRIC_COLUMNS) <= set(metrics.columns)
    assert metrics.output_relpath.dropna().map(lambda v: ':' + chr(92) not in str(v)).all()
