"""Level C2 checks: the hazards the Owner named, and nothing else.

C2-A is a fixed dataset representation, so its guards say what may *not* influence a
sample's ranks: other samples, gene order, or missing values.  C2-B is a model-time
contract, so its guards say what may not be estimated from held-out data.  The last tests
pin the two claims that would be quietest to lose: that no whole-dataset standardized
matrix exists, and that C2 never rates a pair ready which C1 did not.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))

import build_representation_layer as layer
from representation_transforms import (
    RobustFitConfig, apply_robust_standardizer, build_within_sample_rank_matrix,
    fit_robust_standardizer, midrank_percentile, robust_fit_summary)

LAYER = layer.LAYER

GENES = ['G1', 'G2', 'G3', 'G4', 'G5', 'G6']


def sample_frame(values: dict[str, list[float]]) -> pd.DataFrame:
    """One column per sample, one row per gene."""
    return pd.DataFrame(values, index=GENES)


def gene_frame(rows: dict[str, list[float]]) -> pd.DataFrame:
    """Same layout, written the way a C2-B fit reads it: one row of samples per gene."""
    samples = list(next(iter(rows.values())))
    names = [f'sample_{index}' for index in range(len(samples))]
    columns = {name: [rows[gene][index] for gene in GENES] for name, index in
               zip(names, range(len(samples)))}
    return pd.DataFrame(columns, index=GENES)


# --------------------------------------------------------------------------- C2-A


def test_ranking_one_sample_does_not_see_another_sample():
    """Percentiles are within-sample, so no cohort statistic can be hiding in C2-A."""
    base = sample_frame({'only': [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]})
    alone, _ = build_within_sample_rank_matrix(base, GENES)
    with_other, _ = build_within_sample_rank_matrix(
        pd.concat([base, sample_frame({'other': [6.0, 5.0, 4.0, 3.0, 2.0, 1.0]})], axis=1),
        GENES)

    pd.testing.assert_series_equal(alone['only'], with_other['only'], check_dtype=False)


def test_tied_genes_get_one_percentile_and_gene_order_never_breaks_the_tie():
    """A pseudobulk sample can have thousands of zero-valued genes; their ranks must not
    depend on which of them happens to be listed first."""
    tied = sample_frame({'a': [0.0, 0.0, 0.0, 0.0, 8.0, 9.0]})
    ranks, qc = build_within_sample_rank_matrix(tied, GENES)
    block = ranks['a'].loc[['G1', 'G2', 'G3', 'G4']]

    assert block.nunique() == 1
    assert float(block.iloc[0]) == pytest.approx((2.5 - 0.5) / 6.0)
    assert int(qc.finite_genes.iloc[0]) == 6
    assert int(qc.rank_unique_values.iloc[0]) == 3

    # the same numbers, listed the other way round: the tie block is now at the end
    reversed_matrix = tied.loc[list(reversed(GENES))]
    reversed_ranks, _ = build_within_sample_rank_matrix(reversed_matrix, list(reversed(GENES)))
    assert float(reversed_ranks.sort_index().loc['G1', 'a']) == pytest.approx(
        float(block.iloc[0]))
    pd.testing.assert_series_equal(ranks.sort_index()['a'], reversed_ranks.sort_index()['a'],
                                   check_dtype=False)


def test_missing_genes_stay_missing_and_leave_the_denominator():
    """No imputation: NaN genes are neither ranked nor filled with zero, and they do not
    dilute the percentile grid of the genes that are present."""
    values = sample_frame({'a': [1.0, np.nan, 3.0, np.nan, 5.0, 6.0]})
    ranks, qc = build_within_sample_rank_matrix(values, GENES)

    assert np.isnan(ranks.loc[['G2', 'G4'], 'a']).all()
    assert int(qc.finite_genes.iloc[0]) == 4
    assert float(qc.finite_fraction.iloc[0]) == pytest.approx(4 / 6)
    # the denominator is the 4 finite genes, not the 6-gene universe
    assert float(ranks.loc['G1', 'a']) == pytest.approx(0.5 / 4)
    assert float(ranks.loc['G6', 'a']) == pytest.approx(3.5 / 4)


def test_a_strictly_monotone_transform_cannot_change_a_rank():
    """C1 leaves each source on its own scale, so C2-A must be indifferent to the monotone
    difference between those scales, which is exactly what a rank representation is for."""
    linear = sample_frame({'a': [0.0, 1.0, 3.0, 7.0, 15.0, 31.0]})
    plain, _ = build_within_sample_rank_matrix(linear, GENES)
    logged, _ = build_within_sample_rank_matrix(np.log2(linear + 1.0), GENES)

    pd.testing.assert_series_equal(plain['a'], logged['a'], check_dtype=False)


def test_the_published_rank_matrices_use_the_frozen_core_and_nothing_else():
    """The ranking universe is the strict canonical core read from the gene-space artifact:
    one shared coordinate system, not each source's own gene list."""
    core = layer.core_ranking_universe()
    assert len(core) > 1000  # read, never hard-coded, but an empty universe is a failure
    metrics = pd.read_csv(LAYER / 'SOURCE_RANK_METRICS.csv')
    assert metrics.ranking_universe_genes.unique().tolist() == [len(core)]

    qc = pd.read_csv(LAYER / 'SAMPLE_RANK_QC.csv.gz', low_memory=False)
    assert qc.ranking_universe_genes.unique().tolist() == [len(core)]
    assert (qc.finite_genes <= len(core)).all()
    assert qc.finite_fraction.between(0.0, 1.0).all()

    for record in metrics.itertuples():
        assert record.rank_matrix_relpath, record.source_expression
        path = layer.bg.PROJECT / record.rank_matrix_relpath
        if not path.exists():
            continue  # matrices are local-only; the published metric still names its file
        ranks = pd.read_parquet(path)
        assert ranks.index.name == 'hgnc_symbol'
        assert ranks.index.tolist() == core
        assert str(ranks.dtypes.unique()[0]) == 'float32'
        values = ranks.to_numpy(dtype=float)
        values = values[np.isfinite(values)]
        assert (values > 0.0).all() and (values < 1.0).all()


def test_the_rank_formula_is_the_midrank_percentile_the_owner_specified():
    """Strictly inside (0, 1): the extremes are 0.5/N and 1 - 0.5/N, never 0 and 1."""
    ranks = midrank_percentile(pd.Series([3.0, 1.0, 2.0, 100.0, -7.0]))
    assert float(ranks.min()) == pytest.approx(0.5 / 5)
    assert float(ranks.max()) == pytest.approx(1.0 - 0.5 / 5)
    assert list(ranks.round(6)) == [0.7, 0.3, 0.5, 0.9, 0.1]
    assert float(midrank_percentile(pd.Series([4.0])).iloc[0]) == 0.5
    assert np.isnan(midrank_percentile(pd.Series([np.nan])).iloc[0])
    partial = midrank_percentile(pd.Series([np.nan, 2.0, np.nan, 4.0]))
    assert int(partial.notna().sum()) == 2
    assert float(partial.iloc[3]) == pytest.approx(0.75)


# --------------------------------------------------------------------------- C2-B


def test_heldout_values_cannot_move_a_train_fit():
    """The whole point of C2-B: median and MAD see the training samples and nothing else."""
    train = gene_frame({'G1': [1.0, 2.0, 3.0, 4.0], 'G2': [2.0, 3.0, 4.0, 5.0],
                        'G3': [3.0, 4.0, 5.0, 6.0], 'G4': [4.0, 5.0, 6.0, 7.0],
                        'G5': [1.0, 4.0, 7.0, 10.0], 'G6': [2.0, 2.0, 9.0, 9.0]})
    state = fit_robust_standardizer(train, ['sample_0', 'sample_1', 'sample_2', 'sample_3'],
                                    source_expression='TEST', split_id='split_a',
                                    upstream_fold_isolation='FOLD_INDEPENDENT')

    contaminated = train.copy()
    contaminated['sample_4'] = [900.0 + index for index in range(len(GENES))]
    contaminated['sample_5'] = [-900.0 - index for index in range(len(GENES))]
    again = fit_robust_standardizer(contaminated,
                                    ['sample_0', 'sample_1', 'sample_2', 'sample_3'],
                                    source_expression='TEST', split_id='split_a',
                                    upstream_fold_isolation='FOLD_INDEPENDENT')

    pd.testing.assert_frame_equal(state, again)
    assert robust_fit_summary(state)['genes_ready'] == len(GENES)


def test_transform_estimates_nothing_from_what_it_transforms():
    """The same fitted state must give the same numbers whether it is applied to the test
    columns alone or inside the full matrix."""
    train = gene_frame({'G1': [1.0, 2.0, 3.0, 4.0], 'G2': [2.0, 3.0, 4.0, 5.0],
                        'G3': [3.0, 4.0, 5.0, 6.0], 'G4': [4.0, 5.0, 6.0, 7.0],
                        'G5': [1.0, 4.0, 7.0, 10.0], 'G6': [2.0, 2.0, 9.0, 9.0]})
    train_columns = ['sample_0', 'sample_1', 'sample_2', 'sample_3']
    state = fit_robust_standardizer(train, train_columns, source_expression='TEST',
                                    split_id='split_a',
                                    upstream_fold_isolation='FOLD_INDEPENDENT')

    held_out = train[['sample_0']].mul(3.0)
    held_out.columns = ['held_out']
    alone = apply_robust_standardizer(held_out, state, expected_source_expression='TEST',
                                      expected_split_id='split_a')
    together = apply_robust_standardizer(pd.concat([train, held_out], axis=1), state,
                                         expected_source_expression='TEST',
                                         expected_split_id='split_a')

    pd.testing.assert_series_equal(alone['held_out'], together['held_out'], check_dtype=False)
    reference = state.set_index('hgnc_symbol')
    expected = ((held_out['held_out'].to_numpy(dtype=float)
                 - reference.median_train.to_numpy(dtype=float))
                / reference.robust_scale.to_numpy(dtype=float))
    np.testing.assert_allclose(together['held_out'].to_numpy(dtype=float), expected,
                               rtol=1e-5, atol=1e-6)


def test_zero_mad_with_spread_falls_back_to_the_iqr():
    """0 0 0 0 4 7 has no median absolute deviation but is not constant; the IQR is the
    second estimate, not an epsilon."""
    values = gene_frame({'G1': [0.0, 0.0, 0.0, 0.0, 4.0, 7.0],
                         'G2': [0.0, 0.0, 0.0, 1.0, 5.0, 8.0],
                         'G3': [1.0, 1.0, 1.0, 2.0, 3.0, 9.0],
                         'G4': [0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
                         'G5': [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
                         'G6': [2.0, 3.0, 4.0, 5.0, 6.0, 7.0]})
    state = fit_robust_standardizer(values, list(values.columns), source_expression='TEST',
                                    split_id='split_b',
                                    upstream_fold_isolation='FOLD_INDEPENDENT')
    first = state.set_index('hgnc_symbol').loc['G1']

    assert float(first.mad_train) == 0.0
    assert float(first.iqr_train) > 0.0
    assert first.scale_method == 'IQR_FALLBACK'
    assert first.fit_status == 'READY'
    assert float(first.robust_scale) == pytest.approx(
        float(first.iqr_train) / 1.3489795003921634)
    assert state.set_index('hgnc_symbol').loc['G5'].scale_method == 'MAD'


def test_a_gene_with_no_spread_at_all_is_unusable_not_rescaled_by_an_epsilon():
    """Flooring a zero scale at 1e-8 would turn a constant gene into an enormous z-score."""
    values = gene_frame({'G1': [5.0, 5.0, 5.0, 5.0], 'G2': [1.0, 2.0, 3.0, 4.0],
                         'G3': [1.0, 4.0, 7.0, 10.0], 'G4': [0.0, 1.0, 8.0, 9.0],
                         'G5': [3.0, 4.0, 5.0, 9.0], 'G6': [1.0, 2.0, 2.0, 3.0]})
    state = fit_robust_standardizer(values, list(values.columns), source_expression='TEST',
                                    split_id='split_c',
                                    upstream_fold_isolation='FOLD_INDEPENDENT')
    constant = state.set_index('hgnc_symbol').loc['G1']

    assert constant.fit_status == 'UNUSABLE_CONSTANT_OR_SPARSE'
    assert constant.scale_method == 'UNUSABLE'
    assert np.isnan(float(constant.robust_scale))

    output = apply_robust_standardizer(values, state, expected_source_expression='TEST',
                                       expected_split_id='split_c')
    assert np.isnan(output.loc['G1'].to_numpy(dtype=float)).all()
    usable = output.drop(index='G1').to_numpy(dtype=float)
    assert np.isfinite(usable).all()
    # an epsilon floor would have made these |z| enormous instead of ordinary
    assert float(np.abs(usable).max()) < 1e4
    summary = robust_fit_summary(state)
    assert summary['genes_constant_or_sparse'] == 1
    assert summary['genes_ready'] == len(GENES) - 1


def test_sparse_genes_are_not_fitted_from_an_empty_train_slice():
    values = gene_frame({'G1': [1.0, 2.0, 3.0, 4.0], 'G2': [1.0, 2.0, 3.0, 4.0],
                         'G3': [1.0, 2.0, 3.0, 4.0], 'G4': [1.0, 2.0, 3.0, 4.0],
                         'G5': [1.0, 2.0, 3.0, 4.0], 'G6': [1.0, 2.0, 3.0, 4.0]})
    values.loc['G1', ['sample_0', 'sample_1']] = np.nan
    state = fit_robust_standardizer(values, list(values.columns), source_expression='TEST',
                                    split_id='split_d',
                                    upstream_fold_isolation='FOLD_INDEPENDENT',
                                    config=RobustFitConfig(min_train_samples_per_gene=4))

    assert (state.set_index('hgnc_symbol').loc['G1'].fit_status
            == 'INSUFFICIENT_TRAIN_SAMPLES')
    assert int(state.set_index('hgnc_symbol').loc['G1'].n_train_finite) == 2
    assert robust_fit_summary(state)['genes_insufficient_train_samples'] == 1


def test_a_nonisolated_upstream_representation_fails_closed():
    """GSE319641 is AUTHOR_FULL_SOURCE ComBat: its upstream saw every cohort including the
    test patients, so C2-B may not present it as leakage-isolated.  Ranking is
    fold-independent; that source's representation is not."""
    values = six_sample_frame()
    train_columns = list(values.columns)
    with pytest.raises(ValueError, match='refuses to present this as leakage-isolated'):
        fit_robust_standardizer(values, train_columns, source_expression='GSE319641_bulk',
                                split_id='split_e', upstream_fold_isolation='NOT_ESTABLISHED')

    stressed = fit_robust_standardizer(values, train_columns, source_expression='GSE319641_bulk',
                                       split_id='split_e',
                                       upstream_fold_isolation='NOT_ESTABLISHED',
                                       allow_nonisolated_upstream=True)
    assert len(stressed) == len(GENES)

    for column in ([], ['absent_sample'], ['sample_0', 'sample_0']):
        with pytest.raises(ValueError):
            fit_robust_standardizer(values, column, source_expression='TEST', split_id='split_e',
                                    upstream_fold_isolation='FOLD_INDEPENDENT')


def test_forgetting_the_fold_isolation_declaration_does_not_borrow_a_safe_one():
    """The safe default was the hole: a caller who omits the argument used to be read as
    declaring the upstream fold-independent, which is exactly how GSE319641 would have
    slipped past the guard it exists for."""
    values = six_sample_frame()
    train_columns = list(values.columns)
    with pytest.raises(ValueError, match='must be explicitly declared'):
        fit_robust_standardizer(values, train_columns, source_expression='GSE319641_bulk',
                                split_id='fold1')

    with pytest.raises(ValueError, match='must be explicitly declared'):
        fit_robust_standardizer(values, train_columns, source_expression='TEST',
                                split_id='fold1', upstream_fold_isolation=None)

    # declaring it is what makes the call legal, and only for a fold-independent upstream
    declared = fit_robust_standardizer(values, train_columns, source_expression='TEST',
                                       split_id='fold1',
                                       upstream_fold_isolation='FOLD_INDEPENDENT')
    assert len(declared) == len(GENES)


def test_a_scaler_cannot_be_applied_to_another_source_or_another_split():
    """GSE91061's fold-1 median and MAD fitted onto MORRISON is a different transform, and
    an automated pipeline would never notice."""
    values = six_sample_frame()
    train_columns = list(values.columns)
    state = fit_robust_standardizer(values, train_columns, source_expression='GSE91061_bulk',
                                    split_id='fold_1',
                                    upstream_fold_isolation='FOLD_INDEPENDENT')

    with pytest.raises(ValueError, match='was fitted on source_expression=GSE91061_bulk'):
        apply_robust_standardizer(values, state, expected_source_expression='MORRISON_bulk',
                                  expected_split_id='fold_1')
    with pytest.raises(ValueError, match='was fitted on split_id=fold_1'):
        apply_robust_standardizer(values, state, expected_source_expression='GSE91061_bulk',
                                  expected_split_id='fold_2')


def test_a_mixed_source_or_mixed_split_state_is_not_a_scaler():
    values = six_sample_frame()
    train_columns = list(values.columns)
    first = fit_robust_standardizer(values, train_columns, source_expression='SOURCE_A',
                                    split_id='fold_1',
                                    upstream_fold_isolation='FOLD_INDEPENDENT')
    other_source = fit_robust_standardizer(values, train_columns, source_expression='SOURCE_B',
                                           split_id='fold_1',
                                           upstream_fold_isolation='FOLD_INDEPENDENT')
    other_split = fit_robust_standardizer(values, train_columns, source_expression='SOURCE_A',
                                          split_id='fold_2',
                                          upstream_fold_isolation='FOLD_INDEPENDENT')

    for blended, phrase in ((pd.concat([first, other_source]), 'mixes 2 values of '
                                                            'source_expression'),
                            (pd.concat([first, other_split]),
                            'mixes 2 values of split_id')):
        with pytest.raises(ValueError, match=phrase):
            apply_robust_standardizer(values, blended, expected_source_expression='SOURCE_A',
                                      expected_split_id='fold_1')

    anonymous = first.drop(columns=['source_expression'])
    with pytest.raises(ValueError, match='missing columns'):
        apply_robust_standardizer(values, anonymous, expected_source_expression='SOURCE_A',
                                  expected_split_id='fold_1')


def test_a_matrix_with_duplicate_rows_or_columns_is_rejected_before_any_transform():
    """A repeated sample column would be weighted twice, and a repeated gene row makes the
    fitted state ambiguous about which row a median belongs to."""
    values = six_sample_frame()
    state = fit_robust_standardizer(values, list(values.columns), source_expression='TEST',
                                    split_id='split_a',
                                    upstream_fold_isolation='FOLD_INDEPENDENT')
    call = dict(expected_source_expression='TEST', expected_split_id='split_a')

    duplicated_genes = pd.concat([values, values.iloc[:1]], axis=0)
    with pytest.raises(ValueError, match='matrix to transform contains duplicate gene rows'):
        apply_robust_standardizer(duplicated_genes, state, **call)

    duplicated_samples = pd.concat([values, values.iloc[:, :1]], axis=1)
    with pytest.raises(ValueError, match='duplicate sample columns'):
        apply_robust_standardizer(duplicated_samples, state, **call)

    duplicated_state = pd.concat([state, state.iloc[:1]])
    with pytest.raises(ValueError, match='state contains duplicate gene rows'):
        apply_robust_standardizer(values, duplicated_state, **call)


# -------------------------------------------------- builder guards, not builder metrics


def six_sample_frame() -> pd.DataFrame:
    return gene_frame({'G1': [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
                       'G2': [2.0, 3.0, 4.0, 5.0, 6.0, 7.0],
                       'G3': [3.0, 4.0, 5.0, 6.0, 7.0, 8.0],
                       'G4': [4.0, 5.0, 6.0, 7.0, 8.0, 9.0],
                       'G5': [1.0, 4.0, 7.0, 10.0, 13.0, 16.0],
                       'G6': [2.0, 2.0, 9.0, 9.0, 4.0, 4.0]})


def fake_c1_matrix(tmp_path: Path) -> str:
    """A two-sample core-only expression matrix, addressed by absolute path so the test
    never touches the real expression_v0.3 or expression_v0.4_rank directories."""
    matrix = pd.DataFrame({'col_a': [1.0, 2.0, 3.0, 4.0], 'col_b': [4.0, 3.0, 2.0, 1.0]},
                          index=GENES[:4])
    path = tmp_path / 'fake_c1.parquet'
    matrix.to_parquet(path)
    return str(path)


FAKE_BINDING = pd.DataFrame({'source_expression': ['FAKE'] * 2,
                             'sample_id': ['SAMPLE_A', 'SAMPLE_B'],
                             'matrix_column': ['col_a', 'col_b']})
FAKE_CONTRACT = pd.Series(dict(modality_class='BULK_RNA', input_scale='TPM',
                               status='QUANTITATIVE_READY',
                               preprocessing_scope='OUR_C1_PER_SAMPLE',
                               fold_isolation='FOLD_INDEPENDENT',
                               author_batch_corrected=False))


def fake_ranker(cell: float):
    """A ranker that puts a fixed value in the top-left cell, to prove the boundary
    invariant is enforced rather than merely recorded."""

    def build(matrix, ranking_genes):
        ranks = pd.DataFrame(0.5, index=list(ranking_genes), columns=matrix.columns,
                             dtype=np.float32)
        ranks.iloc[0, 0] = cell
        qc = pd.DataFrame({'sample_id': list(matrix.columns),
                           'ranking_universe_genes': len(list(ranking_genes)),
                           'finite_genes': len(list(ranking_genes)),
                           'finite_fraction': 1.0, 'rank_min': 0.1, 'rank_max': 0.9,
                           'rank_unique_values': 4})
        return ranks, qc

    return build


def test_the_builder_refuses_to_write_a_matrix_that_breaks_a_rank_invariant(tmp_path,
                                                                            monkeypatch):
    monkeypatch.setattr(layer, 'build_within_sample_rank_matrix', fake_ranker(1.0))
    with pytest.raises(ValueError, match='touches 0 or 1'):
        layer.rank_one_source('FAKE', fake_c1_matrix(tmp_path), FAKE_BINDING, GENES[:4],
                              FAKE_CONTRACT)

    monkeypatch.setattr(layer, 'build_within_sample_rank_matrix', fake_ranker(float('nan')))
    with pytest.raises(ValueError, match='received no rank'):
        layer.rank_one_source('FAKE', fake_c1_matrix(tmp_path), FAKE_BINDING, GENES[:4],
                              FAKE_CONTRACT)


def test_the_builder_refuses_a_ranker_that_lets_other_samples_move_a_sample(tmp_path,
                                                                           monkeypatch):
    """The leak-shaped hazard for a rank representation: a ranker whose output depends on
    how many samples are in the frame must not survive the build."""

    def leaky(matrix, ranking_genes):
        ranks = pd.DataFrame(0.5, index=list(ranking_genes), columns=matrix.columns,
                             dtype=np.float32)
        if len(matrix.columns) > 1:
            ranks.iloc[0, 0] = 0.9
        qc = pd.DataFrame({'sample_id': list(matrix.columns),
                           'ranking_universe_genes': len(list(ranking_genes)),
                           'finite_genes': len(list(ranking_genes)),
                           'finite_fraction': 1.0, 'rank_min': 0.1, 'rank_max': 0.9,
                           'rank_unique_values': 4})
        return ranks, qc

    monkeypatch.setattr(layer, 'build_within_sample_rank_matrix', leaky)
    with pytest.raises(ValueError, match='moved by'):
        layer.rank_one_source('FAKE', fake_c1_matrix(tmp_path), FAKE_BINDING, GENES[:4],
                              FAKE_CONTRACT)

    assert layer.independence_violation(0.0) == ''
    assert layer.independence_violation(layer.INDEPENDENCE_TOLERANCE) == ''
    assert 'moved by' in layer.independence_violation(0.25)
    assert 'touches 0 or 1' in layer.rank_invariant_violation(np.array([0.25, 1.0]), 2)
    assert 'no finite' in layer.rank_invariant_violation(np.array([]), 0)
    assert ('received no rank' in layer.rank_invariant_violation(np.array([np.nan, 0.5]), 2))
    assert layer.rank_invariant_violation(np.array([0.25, 0.75]), 2) == ''


def test_a_partial_build_does_not_report_success():
    """Writing the diagnostics and then exiting 0 is how a partial layer turns green."""
    empty = pd.DataFrame(columns=['source_expression', 'failure'])
    assert layer.build_incompleteness(empty) == ''

    one_failure = pd.DataFrame([{'source_expression': 'MORRISON_bulk',
                                 'failure': 'rank value touches 0 or 1'}])
    message = layer.build_incompleteness(one_failure)
    assert '1 QUANTITATIVE_READY C1 source(s) failed' in message
    assert 'MORRISON_bulk' in message and 'touches 0 or 1' in message


# ------------------------------------------------------- what must not have been built


def test_no_whole_dataset_standardized_matrix_exists():
    """C2-B ships as a contract.  If a pre-standardized matrix ever appeared, this layer's
    leak-prevention claim would be void, so its absence is asserted rather than trusted."""
    report = json.loads((LAYER / 'REPRESENTATION_LAYER_REPORT.json').read_text(encoding='utf-8'))
    assert report['c2b']['global_standardized_matrix_created'] is False
    assert report['c2b']['status'] == 'CONTRACT_ONLY_NOT_MATERIALIZED'

    assert not [path for path in LAYER.rglob('*')
                if 'standard' in path.name.lower() or path.suffix == '.parquet']
    if layer.RANK_MATRIX_DIR.exists():
        assert not [path.name for path in layer.RANK_MATRIX_DIR.iterdir()
                    if 'standard' in path.name.lower()
                    or path.name.endswith('.canonical_log_expression.parquet')]


def test_c2_never_recovers_a_pair_that_c1_did_not_rate_ready():
    """The published coverage artifact must be a subset of C1 readiness, and its counts must
    agree with the per-sample QC it was built from."""
    coverage = pd.read_csv(LAYER / 'PAIR_RANK_COVERAGE.csv.gz', low_memory=False)
    c1 = pd.read_csv(layer.EXPRESSION_LAYER / 'PAIR_EXPRESSION_COVERAGE.csv.gz', low_memory=False)
    c1_ready = c1[c1.quantitative_status.eq('QUANTITATIVE_READY')
                  & (c1.matrix_column_t0 != c1.matrix_column_t1)]

    ready = coverage[coverage.rank_status.eq(layer.RANK_READY)]
    assert set(ready.pair_uid) <= set(c1_ready.pair_uid)
    assert len(coverage) == len(c1)
    unresolved = coverage[coverage.rank_status.eq(layer.RANK_PARENT_NOT_READY)]
    assert set(unresolved.rank_detail) <= set(c1.quantitative_status.unique())

    qc = pd.read_csv(LAYER / 'SAMPLE_RANK_QC.csv.gz', low_memory=False)
    ranked = {(row.source_expression, row.sample_id) for row in qc.itertuples()
              if row.finite_genes > 0}
    for row in ready.itertuples():
        assert (row.source_expression, row.sample_t0) in ranked
        assert (row.source_expression, row.sample_t1) in ranked
        assert row.fold_isolation in ('FOLD_INDEPENDENT', 'NOT_ESTABLISHED')
    assert set(ready[ready.author_batch_corrected.astype(bool)].source_expression.unique()) \
        <= {'GSE319641_bulk'}
    assert qc[qc.author_batch_corrected.astype(bool)].preprocessing_scope.unique().tolist() \
        == ['AUTHOR_FULL_SOURCE']


def test_matrices_are_reachable_only_through_the_manifest_of_this_build():
    """A failed source leaves its previous parquet on disk, so the directory is not the
    truth: the manifest is.  Every hash is checked against the file it names, which is what
    turns a stale matrix into a detectable condition rather than a silent one."""
    metrics = pd.read_csv(LAYER / 'SOURCE_RANK_METRICS.csv')
    assert {'build_version', 'rank_matrix_sha256', 'rank_matrix_relpath'} \
        <= set(metrics.columns)
    assert metrics.build_version.unique().tolist() == [layer.REPRESENTATION_LAYER_VERSION]
    assert not metrics.source_expression.duplicated().any()

    built = metrics[metrics.rank_matrix_size_bytes.gt(0)]
    assert built.rank_matrix_sha256.str.fullmatch(r'[0-9a-f]{64}').all()

    on_disk = [row for row in built.itertuples()
               if (layer.bg.PROJECT / row.rank_matrix_relpath).exists()]
    if not on_disk:
        return  # the published copy carries the manifest, not the matrices
    for row in on_disk:
        assert layer.file_sha256(layer.bg.PROJECT / row.rank_matrix_relpath) \
            == row.rank_matrix_sha256
    referenced = set(metrics.rank_matrix_relpath)
    assert layer.stale_rank_matrices(referenced) == []
    # the detector is real: with no manifest claim, every local matrix is unreferenced
    assert len(layer.stale_rank_matrices(set())) == len(on_disk)
