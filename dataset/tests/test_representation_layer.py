"""Level C2 checks: the hazards the Owner named, and nothing else.

C2-A's claim is that a sample's percentile comes from that sample's own vector, so these
tests attack it as an outsider would: delete the other samples, make them enormous, shuffle
them, and ask whether the target column moved.  It cannot, because rank_one_sample receives
one Series and has no route to a neighbour.

C2-B's claim is that a scaler was estimated from training membership only, so the tests hand
it the full matrix, another cohort's state, a split with one sample on both sides, and a
gene set that is the right size but not the right genes.

The last tests pin what must be absent: a whole-dataset standardized matrix, a ready pair C1
did not rate, and any route to a stale matrix that does not name the build that made it.
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
import representation_transforms as rt
from representation_transforms import (
    RobustFitConfig, RobustSplit, RobustState, apply_robust_standardizer,
    build_within_sample_rank_matrix, fit_robust_standardizer, rank_one_sample,
    robust_fit_summary, verify_fit_membership)

LAYER = layer.LAYER

GENES = ['G1', 'G2', 'G3', 'G4', 'G5', 'G6']


def sample_frame(values: dict[str, list[float]]) -> pd.DataFrame:
    """One column per sample, one row per gene — the layout C2-A reads."""
    return pd.DataFrame(values, index=GENES)


def six_sample_frame() -> pd.DataFrame:
    return pd.DataFrame(
        {f'sample_{index}': [row[index] for row in (
            (1.0, 2.0, 3.0, 4.0, 5.0, 6.0), (2.0, 3.0, 4.0, 5.0, 6.0, 7.0),
            (3.0, 4.0, 5.0, 6.0, 7.0, 8.0), (4.0, 5.0, 6.0, 7.0, 8.0, 9.0),
            (1.0, 4.0, 7.0, 10.0, 13.0, 16.0), (2.0, 2.0, 9.0, 9.0, 4.0, 4.0))]
         for index in range(6)}, index=GENES)


def fold_split(values: pd.DataFrame, source: str = 'TEST', split_id: str = 'fold_1',
               train: int = 4) -> RobustSplit:
    """A split whose membership is real: the first `train` samples train, the rest evaluate."""
    return RobustSplit(source_expression=source, split_id=split_id,
                       train_sample_ids=list(values.columns[:train]),
                       eval_sample_ids=list(values.columns[train:]),
                       upstream_fold_isolation='FOLD_INDEPENDENT')


def fit_on_train(values: pd.DataFrame, split: RobustSplit, **kwargs) -> RobustState:
    return fit_robust_standardizer(values[list(split.train_sample_ids)], split, **kwargs)


# --------------------------------------------------------------------------- C2-A


def test_deleting_other_samples_cannot_move_a_sample():
    base = sample_frame({'target': [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]})
    alone, _ = build_within_sample_rank_matrix(base, GENES)
    together, _ = build_within_sample_rank_matrix(
        pd.concat([base, sample_frame({'other': [6.0, 5.0, 4.0, 3.0, 2.0, 1.0]})], axis=1),
        GENES)

    pd.testing.assert_series_equal(alone['target'], together['target'], check_dtype=False)


def test_extremizing_or_shuffling_other_samples_cannot_move_a_sample():
    """The two ways a cohort statistic hides inside a per-sample transform: a column that
    dominates every gene, and a different column order."""
    target = sample_frame({'target': [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]})
    neighbours = {'a': [1.0] * len(GENES), 'b': [2.0] * len(GENES)}
    ordinary, _ = build_within_sample_rank_matrix(
        pd.concat([target, sample_frame(neighbours)], axis=1), GENES)
    enormous, _ = build_within_sample_rank_matrix(
        pd.concat([target, sample_frame({'a': [1e9] * len(GENES),
                                         'b': [-1e9] * len(GENES)})], axis=1), GENES)
    shuffled, _ = build_within_sample_rank_matrix(
        pd.concat([target, sample_frame(neighbours)], axis=1)[['b', 'a', 'target']], GENES)

    pd.testing.assert_series_equal(ordinary['target'], enormous['target'], check_dtype=False)
    pd.testing.assert_series_equal(ordinary['target'], shuffled['target'], check_dtype=False)
    assert float(enormous['target'].min()) == pytest.approx(0.5 / len(GENES))
    assert float(enormous['target'].max()) == pytest.approx(1.0 - 0.5 / len(GENES))


def test_the_ranking_primitive_receives_one_sample_and_nothing_else():
    """The spy proves the loop that assembles a matrix never hands that primitive more than
    one sample, and the equality proves nothing else contributes to a column either."""
    values = sample_frame({'a': [3.0, 1.0, 2.0, 8.0, 5.0, 4.0],
                           'b': [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
                           'c': [6.0, 5.0, 4.0, 3.0, 2.0, 1.0]})

    seen: list = []
    real = rt.rank_one_sample

    def spy(gene_values):
        seen.append(gene_values)
        return real(gene_values)

    rt.rank_one_sample = spy
    try:
        matrix, _ = build_within_sample_rank_matrix(values, GENES)
    finally:
        rt.rank_one_sample = real

    assert len(seen) == values.shape[1]
    for argument in seen:
        assert isinstance(argument, pd.Series), 'the primitive was handed a whole matrix'
        assert len(argument) == len(GENES)
        assert set(argument.index) == set(GENES)
    for sample in values.columns:
        pd.testing.assert_series_equal(rank_one_sample(values[sample]), matrix[sample],
                                       check_dtype=False, check_names=False)


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

    reversed_ranks, _ = build_within_sample_rank_matrix(tied.loc[list(reversed(GENES))],
                                                       list(reversed(GENES)))
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
    assert float(ranks.loc['G1', 'a']) == pytest.approx(0.5 / 4)
    assert float(ranks.loc['G6', 'a']) == pytest.approx(3.5 / 4)


def test_a_strictly_monotone_transform_cannot_change_a_rank():
    """C1 leaves each source on its own scale, so C2-A must be indifferent to the monotone
    difference between those scales, which is exactly what a rank representation is for."""
    linear = sample_frame({'a': [0.0, 1.0, 3.0, 7.0, 15.0, 31.0]})
    plain, _ = build_within_sample_rank_matrix(linear, GENES)
    logged, _ = build_within_sample_rank_matrix(np.log2(linear + 1.0), GENES)

    pd.testing.assert_series_equal(plain['a'], logged['a'], check_dtype=False)


def test_the_rank_formula_is_the_midrank_percentile_the_owner_specified():
    """Strictly inside (0, 1): the extremes are 0.5/N and 1 - 0.5/N, never 0 and 1."""
    ranks = rank_one_sample(pd.Series([3.0, 1.0, 2.0, 100.0, -7.0]))
    assert float(ranks.min()) == pytest.approx(0.5 / 5)
    assert float(ranks.max()) == pytest.approx(1.0 - 0.5 / 5)
    assert list(ranks.round(6)) == [0.7, 0.3, 0.5, 0.9, 0.1]
    assert float(rank_one_sample(pd.Series([4.0])).iloc[0]) == 0.5
    assert np.isnan(rank_one_sample(pd.Series([np.nan])).iloc[0])
    partial = rank_one_sample(pd.Series([np.nan, 2.0, np.nan, 4.0]))
    assert int(partial.notna().sum()) == 2
    assert float(partial.iloc[3]) == pytest.approx(0.75)


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


# --------------------------------------------------------------------------- C2-B


def test_a_fit_handed_the_full_matrix_is_refused():
    """The bad implementation this guard exists for: fit over every column and call the
    training fold a subset.  The resulting state would still look like a train fit."""
    values = six_sample_frame()
    split = fold_split(values)

    with pytest.raises(ValueError, match='may only be handed its own training samples'):
        fit_robust_standardizer(values, split)

    state = fit_on_train(values, split)
    assert state.fit_sample_ids == tuple(split.train_sample_ids)
    assert verify_fit_membership(state, split) == ''

    scaled = values.copy()
    for column in split.eval_sample_ids:
        scaled[column] = scaled[column] * 1000.0 + 500.0
    pd.testing.assert_frame_equal(state.frame, fit_on_train(scaled, split).frame)


def test_a_split_cannot_put_one_sample_on_both_sides():
    """Disjointness is constructed here rather than asserted later, which is what makes the
    membership comparison in verify_fit_membership a complete leakage proof."""
    with pytest.raises(ValueError, match='both training and evaluation'):
        RobustSplit(source_expression='TEST', split_id='fold_1',
                    train_sample_ids=['sample_0', 'sample_1'],
                    eval_sample_ids=['sample_1', 'sample_2'])
    with pytest.raises(ValueError, match='no training sample'):
        RobustSplit(source_expression='TEST', split_id='fold_1', train_sample_ids=[])
    with pytest.raises(ValueError, match='train sample ids contain duplicates'):
        RobustSplit(source_expression='TEST', split_id='fold_1',
                    train_sample_ids=['sample_0', 'sample_0'])
    with pytest.raises(ValueError, match='name both its source and its split'):
        RobustSplit(source_expression='', split_id='fold_1', train_sample_ids=['sample_0'])

    values = six_sample_frame()
    split = fold_split(values)
    assert split.known_sample_ids == set(values.columns)


def test_a_scaler_cannot_be_applied_to_another_source_or_another_split():
    """GSE91061's fold-1 median and MAD applied to MORRISON is a different transform, and
    an automated pipeline would never notice."""
    values = six_sample_frame()
    state = fit_on_train(values, fold_split(values, 'GSE91061_bulk', 'fold_1'))
    train = values[list(state.fit_sample_ids)]

    for wrong, phrase in ((fold_split(values, 'MORRISON_bulk', 'fold_1'),
                           'belongs to source GSE91061_bulk'),
                          (fold_split(values, 'GSE91061_bulk', 'fold_2'),
                           'belongs to split fold_1')):
        with pytest.raises(ValueError, match=phrase):
            apply_robust_standardizer(train, state, wrong)

    smuggled = RobustState(frame=state.frame, source_expression='GSE91061_bulk',
                           split_id='fold_1',
                           fit_sample_ids=state.fit_sample_ids + ('sample_4',),
                           fit_gene_ids=state.fit_gene_ids)
    with pytest.raises(ValueError, match='not fitted on exactly this split'):
        apply_robust_standardizer(train, smuggled, fold_split(values, 'GSE91061_bulk',
                                                              'fold_1'))


def test_an_unknown_or_nonisolated_upstream_fails_closed():
    """A split that did not load its fold isolation knows nothing, and silence is not the
    safe answer: only FOLD_INDEPENDENT, carried from the C1 manifest, lets a fit proceed."""
    values = six_sample_frame()
    train = values.iloc[:, :4]

    for isolation in ('', 'NOT_ESTABLISHED', 'AUTHOR_FULL_SOURCE', 'FOLD-FREE'):
        split = RobustSplit(source_expression='TEST', split_id='fold_1',
                            train_sample_ids=list(values.columns[:4]),
                            eval_sample_ids=list(values.columns[4:]),
                            upstream_fold_isolation=isolation)
        with pytest.raises(ValueError, match='refuses to present this as leakage-isolated'):
            fit_robust_standardizer(train, split)

    assert fit_robust_standardizer(
        train, RobustSplit(source_expression='TEST', split_id='fold_1',
                           train_sample_ids=list(values.columns[:4]),
                           eval_sample_ids=list(values.columns[4:]),
                           upstream_fold_isolation='FOLD_INDEPENDENT')).fit_gene_ids == tuple(GENES)

    stressed = fit_robust_standardizer(train, RobustSplit(
        source_expression='GSE319641_bulk', split_id='fold_1',
        train_sample_ids=list(values.columns[:4]), eval_sample_ids=list(values.columns[4:]),
        upstream_fold_isolation='NOT_ESTABLISHED'), allow_nonisolated_upstream=True)
    assert len(stressed.frame) == len(GENES)


def test_gse319641s_isolation_comes_from_the_c1_manifest_not_from_the_caller():
    """The one argument that decides whether a scaler may claim to be leakage-isolated is
    read off the C1 artifact by robust_split_for_source, so a caller cannot type the
    convenient value.  That source's C1 upstream is AUTHOR_FULL_SOURCE ComBat: the authors'
    batch correction already saw the test patients."""
    metrics = pd.read_csv(layer.EXPRESSION_LAYER / 'SOURCE_EXPRESSION_METRICS.csv')
    nonisolated = metrics.loc[~metrics.fold_isolation.astype(str).eq('FOLD_INDEPENDENT'),
                              'source_expression'].tolist()
    binding = pd.read_csv(layer.EXPRESSION_LAYER / 'manifests/SAMPLE_COLUMN_BINDING.csv.gz',
                          low_memory=False)

    def samples_of(source):
        return binding.loc[binding.source_expression.eq(source)
                           & binding.binding_status.eq('BOUND'), 'sample_id'].tolist()

    for source in nonisolated:
        row = metrics[metrics.source_expression.eq(source)]
        assert row.preprocessing_scope.iloc[0] == 'AUTHOR_FULL_SOURCE'
        assert row.status.iloc[0] == 'QUANTITATIVE_READY'
        samples = samples_of(source)
        if len(samples) < 4:
            continue
        split = layer.robust_split_for_source(source, 'fold_1', samples[:2], samples[2:4])
        assert split.upstream_fold_isolation == 'NOT_ESTABLISHED'
        with pytest.raises(ValueError, match='refuses to present this as leakage-isolated'):
            fit_on_train(pd.DataFrame(index=GENES, columns=samples[:2], dtype=float), split)

    isolated = layer.robust_split_for_source('GSE91061_bulk', 'fold_1',
                                             samples_of('GSE91061_bulk')[:4])
    assert isolated.upstream_fold_isolation == 'FOLD_INDEPENDENT'

    with pytest.raises(ValueError, match='not a C1 source'):
        layer.robust_split_for_source('NOT_A_SOURCE', 'fold_1', ['a', 'b'])
    not_ready = metrics.loc[~metrics.status.eq('QUANTITATIVE_READY'), 'source_expression']
    if len(not_ready):
        with pytest.raises(ValueError, match='is not QUANTITATIVE_READY in C1'):
            layer.robust_split_for_source(str(not_ready.iloc[0]), 'fold_1', ['a', 'b'])


def test_transform_estimates_nothing_from_what_it_transforms():
    """The same fitted state must give the same numbers whether it is applied to one test
    column alone or inside the whole matrix this split accounts for."""
    values = six_sample_frame()
    split = fold_split(values)
    state = fit_on_train(values, split)

    held_out = values[['sample_4']].mul(3.0)
    held_out.columns = ['held_out']
    evaluation_only = RobustSplit(source_expression='TEST', split_id='fold_1',
                                  train_sample_ids=split.train_sample_ids,
                                  eval_sample_ids=('held_out',))
    alone = apply_robust_standardizer(held_out, state, evaluation_only)
    together = apply_robust_standardizer(
        pd.concat([values[list(split.train_sample_ids)], held_out], axis=1), state,
        evaluation_only)

    pd.testing.assert_series_equal(alone['held_out'], together['held_out'], check_dtype=False)
    reference = state.frame.set_index('hgnc_symbol')
    expected = ((held_out['held_out'].to_numpy(dtype=float)
                 - reference.median_train.to_numpy(dtype=float))
                / reference.robust_scale.to_numpy(dtype=float))
    np.testing.assert_allclose(together['held_out'].to_numpy(dtype=float), expected,
                               rtol=1e-5, atol=1e-6)


def test_membership_is_compared_by_identity_so_swapping_ids_is_caught():
    """The corruption a count check cannot see: one gene deleted and another added, the
    matrix the same shape, the state fitted on a different measurement."""
    values = six_sample_frame()
    split = fold_split(values)
    state = fit_on_train(values, split)
    train = values[list(split.train_sample_ids)]

    assert layer.membership_check('ids', ['a', 'b'], ['a', 'b']) == ''
    assert 'missing=1[b]' in layer.membership_check('ids', ['a', 'b'], ['a', 'c'])
    assert 'unexpected=1[c]' in layer.membership_check('ids', ['a', 'b'], ['a', 'c'])
    assert 'duplicated=1[a]' in layer.membership_check('ids', ['a', 'b'], ['a', 'a', 'b'])

    for broken, phrase in ((train.rename(index={'G3': 'G7'}),
                            'genes to transform are not the genes that were fitted'),
                           (train.drop(index='G3'),
                            'genes to transform are not the genes that were fitted'),
                           (pd.concat([train, train.iloc[:, :1]], axis=1),
                            'duplicate sample columns')):
        with pytest.raises(ValueError, match=phrase):
            apply_robust_standardizer(broken, state, split)

    for broken in (train.rename(columns={'sample_0': 'sample_9'}), train.drop(columns='sample_0')):
        with pytest.raises(ValueError, match='not exactly that membership'):
            fit_robust_standardizer(broken, split)

    stranger = pd.concat([train, values[['sample_4']].rename(
        columns={'sample_4': 'another_cohorts_sample'})], axis=1)
    with pytest.raises(ValueError, match='neither this split'):
        apply_robust_standardizer(stranger, state, split)


def test_a_matrix_that_doubles_a_sample_or_a_gene_is_not_transformable():
    """A repeated column would weight one patient twice inside a fold that claims to count
    each of them once."""
    values = six_sample_frame()
    split = fold_split(values)
    train = values[list(split.train_sample_ids)]

    with pytest.raises(ValueError, match='duplicate sample columns'):
        fit_robust_standardizer(pd.concat([train, train.iloc[:, :1]], axis=1), split)
    with pytest.raises(ValueError, match='duplicate gene rows'):
        fit_robust_standardizer(pd.concat([train, train.iloc[:1]], axis=0), split)

    state = fit_robust_standardizer(train, split)
    with pytest.raises(ValueError, match='duplicate sample columns'):
        apply_robust_standardizer(pd.concat([train, train.iloc[:, :1]], axis=1), state, split)
    with pytest.raises(ValueError, match='duplicate gene rows'):
        apply_robust_standardizer(pd.concat([train, train.iloc[:1]], axis=0), state, split)
    with pytest.raises(ValueError, match='state contains duplicate gene rows'):
        apply_robust_standardizer(train, RobustState(
            frame=pd.concat([state.frame, state.frame.iloc[:1]]),
            source_expression='TEST', split_id='fold_1',
            fit_sample_ids=state.fit_sample_ids,
            fit_gene_ids=state.fit_gene_ids + ('G1',)), split)

    with pytest.raises(ValueError, match='no training sample'):
        RobustSplit(source_expression='TEST', split_id='fold_1', train_sample_ids=[],
                    upstream_fold_isolation='FOLD_INDEPENDENT')


def test_zero_mad_with_spread_falls_back_to_the_iqr():
    """A training fold of 0 0 0 4 has no median absolute deviation but is not constant; the
    IQR is the second estimate, not an epsilon."""
    values = pd.DataFrame(
        {f'sample_{index}': [row[index] for row in (
            (0.0, 0.0, 0.0, 4.0, 7.0, 9.0), (0.0, 0.0, 1.0, 5.0, 6.0, 8.0),
            (1.0, 1.0, 1.0, 2.0, 3.0, 9.0), (0.0, 1.0, 2.0, 3.0, 4.0, 5.0),
            (1.0, 2.0, 3.0, 4.0, 5.0, 6.0), (2.0, 3.0, 4.0, 5.0, 6.0, 7.0))]
         for index in range(6)}, index=GENES)
    state = fit_on_train(values, fold_split(values))
    frame = state.frame.set_index('hgnc_symbol')
    first = frame.loc['G1']

    assert float(first.mad_train) == 0.0
    assert float(first.iqr_train) > 0.0
    assert first.scale_method == 'IQR_FALLBACK'
    assert first.fit_status == 'READY'
    assert float(first.robust_scale) == pytest.approx(
        float(first.iqr_train) / 1.3489795003921634)
    assert frame.loc['G5'].scale_method == 'MAD'


def test_a_gene_with_no_spread_at_all_is_unusable_not_rescaled_by_an_epsilon():
    """Flooring a zero scale at 1e-8 would turn a constant gene into an enormous z-score."""
    values = six_sample_frame().iloc[:, :4]
    values.loc['G1'] = 5.0
    split = RobustSplit(source_expression='TEST', split_id='fold_1',
                        train_sample_ids=list(values.columns),
                        upstream_fold_isolation='FOLD_INDEPENDENT')
    state = fit_robust_standardizer(values, split)
    constant = state.frame.set_index('hgnc_symbol').loc['G1']

    assert constant.fit_status == 'UNUSABLE_CONSTANT_OR_SPARSE'
    assert constant.scale_method == 'UNUSABLE'
    assert np.isnan(float(constant.robust_scale))

    output = apply_robust_standardizer(values, state, split)
    assert np.isnan(output.loc['G1'].to_numpy(dtype=float)).all()
    usable = output.drop(index='G1').to_numpy(dtype=float)
    assert np.isfinite(usable).all()
    assert float(np.abs(usable).max()) < 1e4
    summary = robust_fit_summary(state)
    assert summary['genes_constant_or_sparse'] == 1
    assert summary['genes_ready'] == len(GENES) - 1


def test_sparse_genes_are_not_fitted_from_an_empty_train_slice():
    values = six_sample_frame()
    values.loc['G1', ['sample_0', 'sample_1', 'sample_2']] = np.nan
    split = RobustSplit(source_expression='TEST', split_id='fold_1',
                        train_sample_ids=list(values.columns),
                        upstream_fold_isolation='FOLD_INDEPENDENT')
    state = fit_robust_standardizer(values, split,
                                    config=RobustFitConfig(min_train_samples_per_gene=4))
    frame = state.frame.set_index('hgnc_symbol')

    assert frame.loc['G1'].fit_status == 'INSUFFICIENT_TRAIN_SAMPLES'
    assert int(frame.loc['G1'].n_train_finite) == 3
    assert np.isnan(float(frame.loc['G1'].median_train))
    assert robust_fit_summary(state)['genes_insufficient_train_samples'] == 1


# -------------------------------------------------- builder guards, not builder metrics


def fake_c1_matrix(directory: Path, genes: list[str]) -> str:
    """A two-sample core-only expression matrix, written under the test's own directory so
    it never touches the real expression_v0.3 matrices or the published rank directory."""
    matrix = pd.DataFrame({'col_a': [1.0, 2.0, 3.0, 4.0], 'col_b': [4.0, 3.0, 2.0, 1.0]},
                          index=genes)
    path = directory / 'fake_c1.parquet'
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
CORE = GENES[:4]
FAKE_RANK_MATRIX = 'FAKE.canonical_core_rank.parquet'


def run_rank_one_source(monkeypatch, tmp_path, ranker, genes=CORE, binding=FAKE_BINDING):
    """Drive the builder over a fake source, with every file it writes landing in tmp_path.

    A manifest locator is project-relative, so the fake project root has to be the test's
    own directory; that also keeps this working in the published copy, which carries the
    manifest but not the matrices.
    """
    output = tmp_path / 'rank'
    output.mkdir(exist_ok=True)
    monkeypatch.setattr(layer, 'build_within_sample_rank_matrix', ranker)
    monkeypatch.setattr(layer, 'RANK_MATRIX_DIR', output)
    monkeypatch.setattr(layer.bg, 'PROJECT', tmp_path)
    return layer.rank_one_source('FAKE', fake_c1_matrix(tmp_path, list(genes)), binding,
                                 list(genes), FAKE_CONTRACT, 'BUILD_TEST')


def fake_ranker(break_it):
    """Give it the honest ranks and then break one property, so the guard under test is the
    builder's check rather than the mathematics."""

    def build(matrix, ranking_genes):
        ranks, qc = rt.build_within_sample_rank_matrix(matrix, ranking_genes)
        return break_it(ranks.copy(), qc)

    return build


def test_the_builder_refuses_to_write_a_matrix_that_breaks_a_rank_invariant(monkeypatch,
                                                                            tmp_path):
    def touch_one(ranks, qc):
        ranks.iloc[0, 0] = 1.0
        return ranks, qc

    with pytest.raises(ValueError, match='touches 0 or 1'):
        run_rank_one_source(monkeypatch, tmp_path, fake_ranker(touch_one))

    def drop_a_value(ranks, qc):
        ranks.iloc[0, 0] = np.nan
        return ranks, qc

    with pytest.raises(ValueError, match='received no rank'):
        run_rank_one_source(monkeypatch, tmp_path, fake_ranker(drop_a_value))

    def all_nan(ranks, qc):
        return ranks * np.nan, qc

    with pytest.raises(ValueError, match='no finite rank value'):
        run_rank_one_source(monkeypatch, tmp_path, fake_ranker(all_nan))

    assert not (tmp_path / 'rank' / FAKE_RANK_MATRIX).exists()
    assert layer.rank_invariant_violation(np.array([0.25, 0.75]), 2) == ''


def test_the_builder_refuses_a_rank_matrix_whose_ids_are_not_the_ones_it_claims(monkeypatch,
                                                                                tmp_path):
    """Right shape, wrong content: a gene row replaced by another gene, a gene row dropped,
    a sample renamed.  A count comparison passes all three."""
    def swap_gene(ranks, qc):
        return ranks.rename(index={ranks.index[0]: 'NOT_IN_THE_CORE'}), qc

    with pytest.raises(ValueError, match='gene rows does not match'):
        run_rank_one_source(monkeypatch, tmp_path, fake_ranker(swap_gene))

    def lose_a_gene(ranks, qc):
        return ranks.drop(ranks.index[0]), qc

    with pytest.raises(ValueError, match='gene rows does not match'):
        run_rank_one_source(monkeypatch, tmp_path, fake_ranker(lose_a_gene))

    def swap_sample(ranks, qc):
        return ranks.rename(columns={'SAMPLE_B': 'SAMPLE_C'}), qc

    with pytest.raises(ValueError, match='sample columns does not match'):
        run_rank_one_source(monkeypatch, tmp_path, fake_ranker(swap_sample))

    def add_a_sample(ranks, qc):
        return pd.concat([ranks, ranks.iloc[:, :1].rename(columns={'SAMPLE_A': 'SAMPLE_C'})],
                         axis=1), qc

    with pytest.raises(ValueError, match='sample columns does not match'):
        run_rank_one_source(monkeypatch, tmp_path, fake_ranker(add_a_sample))

    assert not (tmp_path / 'rank' / FAKE_RANK_MATRIX).exists()

    record, qc = run_rank_one_source(monkeypatch, tmp_path,
                                     rt.build_within_sample_rank_matrix)
    assert record['build_id'] == 'BUILD_TEST'
    assert record['failure'] == ''
    assert qc['sample_id'].tolist() == ['SAMPLE_A', 'SAMPLE_B']
    # only a source that passed its own checks reaches the disk
    assert (tmp_path / 'rank' / FAKE_RANK_MATRIX).exists()


def test_a_healthy_build_publishes_and_a_partial_one_never_does(tmp_path):
    """healthy build A -> source lost on rebuild B -> A's files are still on disk -> a
    consumer asking for B is refused, and is never handed A instead."""
    metrics = pd.read_csv(LAYER / 'SOURCE_RANK_METRICS.csv')
    assert {'build_id', 'rank_matrix_relpath', 'failure'} <= set(metrics.columns)
    # whether the *build* completed is one fact about the build, not 19 facts about its rows
    assert 'build_status' not in metrics.columns
    # governance by hash was deleted on purpose: identity is readable information
    assert not [column for column in metrics.columns
                if 'sha' in column.lower() or 'hash' in column.lower()]
    assert not any(name in dir(layer) for name in ('file_sha256', 'stale_rank_matrices',
                                                   'current_rank_matrices',
                                                   'INDEPENDENCE_TOLERANCE'))

    published = metrics.build_id.unique().tolist()
    assert len(published) == 1
    published = published[0]
    assert len(layer.rank_matrices_for_build(published)) == len(metrics)

    directory = tmp_path / 'representation_layer'
    directory.mkdir()
    metrics.to_csv(directory / 'SOURCE_RANK_METRICS.csv', index=False)

    def report(build_id, status):
        (directory / 'REPRESENTATION_LAYER_REPORT.json').write_text(
            json.dumps({'headline': {'build_id': build_id, 'build_status': status}}),
            encoding='utf-8')

    report(published, layer.BUILD_COMPLETE)
    assert len(layer.rank_matrices_for_build(published, directory)) == len(metrics)

    report('BUILD_B', layer.BUILD_INCOMPLETE)
    for asked in ('BUILD_B', published):
        with pytest.raises(SystemExit, match='never published'):
            layer.rank_matrices_for_build(asked, directory)

    # a later complete build is the published one, and the earlier build's files are not
    # offered to a consumer that named the earlier build
    report('BUILD_B', layer.BUILD_COMPLETE)
    with pytest.raises(SystemExit, match='will not fall back'):
        layer.rank_matrices_for_build(published, directory)
    with pytest.raises(SystemExit, match='no published C2-A build'):
        layer.rank_matrices_for_build('BUILD_A', directory)

    report('BUILD_B', layer.BUILD_COMPLETE)
    metrics[metrics.source_expression != metrics.source_expression.iloc[0]].assign(
        build_id='BUILD_B').to_csv(directory / 'SOURCE_RANK_METRICS.csv', index=False)
    with pytest.raises(SystemExit, match='incomplete on disk'):
        layer.rank_matrices_for_build('BUILD_B', directory)


def test_a_partial_build_does_not_report_success():
    """Writing the diagnostics and then exiting 0 is how a partial layer turns green."""
    assert layer.build_incompleteness(pd.DataFrame(columns=['source_expression',
                                                            'failure'])) == ''
    message = layer.build_incompleteness(pd.DataFrame([
        {'source_expression': 'MORRISON_bulk', 'failure': 'rank value touches 0 or 1'}]))
    assert '1 QUANTITATIVE_READY C1 source(s) failed' in message
    assert 'MORRISON_bulk' in message and 'touches 0 or 1' in message


# ------------------------------------------------------- what must not have been built


def test_no_whole_dataset_standardized_matrix_exists():
    """C2-B ships as a contract.  If a pre-standardized matrix ever appeared, this layer's
    leak-prevention claim would be void, so its absence is asserted rather than trusted."""
    report = json.loads((LAYER / 'REPRESENTATION_LAYER_REPORT.json').read_text(encoding='utf-8'))
    assert report['c2b']['global_standardized_matrix_created'] is False
    assert report['c2b']['status'] == 'CONTRACT_ONLY_NOT_MATERIALIZED'
    assert report['headline']['build_status'] == layer.BUILD_COMPLETE

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
    c1 = pd.read_csv(layer.EXPRESSION_LAYER / 'PAIR_EXPRESSION_COVERAGE.csv.gz',
                     low_memory=False)
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
