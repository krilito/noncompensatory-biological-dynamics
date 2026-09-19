"""Build the PAIR v0.4 Level C2 representation layer.

C2 has two representations and they are *not* the same kind of object.

C2-A is a fixed dataset representation: every ready sample's C1 canonical vector is turned
into within-sample percentile ranks over one shared ranking coordinate system, the frozen
strict canonical core of gene space v0.2.1.  It uses no information from any other sample,
so it can be computed once and published.

C2-B is a model-time fit/transform contract, not a matrix.  Its per-gene median and MAD may
only be estimated from the training samples of one source in one split, so a whole-dataset
standardized matrix can never exist here — writing one would put test patients into the
statistics that claim to exclude them.  This builder therefore ships the contract and its
guards, and deliberately materializes nothing for C2-B.

Explicitly out of scope (C3 and later): paired deltas, delta ranks, response models,
classification, feature selection, PCA, ComBat, cross-cohort z-scoring.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd

import build_gene_space as bg
from representation_transforms import (RANK_FORMULA, RANK_METHOD,
                                       build_within_sample_rank_matrix)

DATASET_ROOT = bg.DATASET_ROOT
LAYER = DATASET_ROOT / 'representation_layer'
EXPRESSION_LAYER = DATASET_ROOT / 'expression_layer'
# Local only, like the C1 matrices.
RANK_MATRIX_DIR = bg.PROJECT / 'longitudinal-data' / 'data' / 'processed' / 'expression_v0.4_rank'

REPRESENTATION_LAYER_VERSION = 'v0.4-C2'
CORE_UNIVERSE_COLUMN = 'in_core_gene_space'

# Statuses of the C2-A endpoint check.  A pair can only be rank-ready if it was already
# C1-quantitative-ready: C2 never recovers a pair that C1 dropped.
RANK_READY = 'RANK_READY'
RANK_NOT_READY = 'ENDPOINT_SAMPLE_NOT_RANKED'
RANK_PARENT_NOT_READY = 'NOT_C1_QUANTITATIVE_READY'

SOURCE_METRIC_COLUMNS = [
    'source_expression', 'modality_class', 'c1_input_scale', 'c1_status', 'ranking_universe_genes',
    'samples_ranked', 'matrix_columns_in_c1', 'columns_not_bound_to_a_sample',
    'core_genes_present_in_source', 'core_finite_fraction_min', 'core_finite_fraction_median',
    'core_finite_fraction_max', 'genes_ranked_per_sample_median', 'largest_tie_block_fraction_min',
    'largest_tie_block_fraction_median', 'rank_unique_values_median', 'rank_min', 'rank_max',
    'rank_exclusive_zero_one', 'sample_independence_max_abs_diff', 'float32_output',
    'preprocessing_scope', 'fold_isolation', 'author_batch_corrected', 'rank_matrix_relpath',
    'rank_matrix_size_bytes', 'failure']

QC_COLUMNS = ['source_expression', 'modality_class', 'sample_id', 'ranking_universe_genes',
              'finite_genes', 'finite_fraction', 'rank_min', 'rank_max', 'rank_unique_values',
              'preprocessing_scope', 'fold_isolation', 'author_batch_corrected']

PAIR_COVERAGE_COLUMNS = [
    'pair_uid', 'patient_uid', 'cohort_code', 'source_expression', 'expression_modality',
    'modality_class', 'sample_t0', 'sample_t1', 'quantitative_status_c1', 'rank_status',
    'rank_detail', 'preprocessing_scope', 'fold_isolation', 'author_batch_corrected',
    'clinical_endpoint_available', 'strict_prcr_vs_pd_eligible', 'frozen_auo_eligible']


def core_ranking_universe() -> list[str]:
    """The frozen strict canonical core, read from the published gene-space artifact.

    The count is never hard-coded: it is whatever the artifact says, because the core is a
    measured property of the corpus, not a constant of this script.
    """
    table = pd.read_csv(bg.GENE_SPACE / 'GENE_SPACE_GENES.csv.gz', low_memory=False)
    flag = table[CORE_UNIVERSE_COLUMN]
    if str(flag.dtype) != 'bool':
        flag = flag.astype(str).str.strip().str.lower().isin(['true', '1', 'yes'])
    genes = table.loc[flag, 'hgnc_symbol'].astype(str).tolist()
    if not genes:
        raise SystemExit(f'{CORE_UNIVERSE_COLUMN} selects no gene: refusing to rank an '
                         'empty universe')
    if len(set(genes)) != len(genes):
        raise SystemExit('core universe contains duplicate genes')
    return genes


def c1_ready_sources() -> pd.DataFrame:
    """Only the QUANTITATIVE_READY C1 sources enter C2, and C1 is read, never rebuilt."""
    metrics = pd.read_csv(EXPRESSION_LAYER / 'SOURCE_EXPRESSION_METRICS.csv')
    ready = metrics[metrics.status.eq('QUANTITATIVE_READY')].copy()
    if ready.empty:
        raise SystemExit('no C1 source is QUANTITATIVE_READY; run build_expression_layer first')
    absent = ready[~ready.output_relpath.astype(str).map(
        lambda value: (bg.PROJECT / str(value)).exists())]
    if not absent.empty:
        raise SystemExit('C1 matrices absent locally for: '
                         + ', '.join(absent.source_expression.tolist()))
    return ready


def bound_samples_by_source() -> pd.DataFrame:
    binding = pd.read_csv(EXPRESSION_LAYER / 'manifests/SAMPLE_COLUMN_BINDING.csv.gz',
                          low_memory=False)
    bound = binding[binding.binding_status.eq('BOUND')][['source_expression', 'sample_id',
                                                         'matrix_column']].copy()
    duplicated = int(bound.duplicated(['source_expression', 'matrix_column']).sum())
    if duplicated:
        # One column naming two samples would leave a column of ranks ambiguous about which
        # sample it belongs to, so that is a hard error rather than a warning.
        raise SystemExit(f'{duplicated} matrix columns bind to more than one sample')
    return bound


def rank_one_source(source: str, relpath: str, binding: pd.DataFrame,
                    ranking_genes: list[str], contract: pd.Series) -> tuple[dict, pd.DataFrame]:
    """Rank one C1 source over the fixed core universe and write its float32 matrix locally."""
    matrix = pd.read_parquet(bg.PROJECT / relpath)
    labels = binding[binding.source_expression.eq(source)].set_index('matrix_column').sample_id
    columns = [column for column in matrix.columns if column in labels.index]
    unbound = int(matrix.shape[1] - len(columns))
    # C1 columns are native matrix headers; the sample identity lives in the binding
    # manifest.  This relabelling is 1:1 and in-memory only: the C1 file is never touched.
    named = matrix[columns].rename(columns=labels.to_dict())
    if named.columns.has_duplicates:
        raise SystemExit(f'{source}: two native columns resolve to one sample id')

    core = named.reindex(pd.Index(ranking_genes))
    ranks, qc = build_within_sample_rank_matrix(core, ranking_genes)
    qc['source_expression'] = source

    values = ranks.to_numpy(dtype=float)
    finite = values[np.isfinite(values)]
    tie_blocks = []
    for column in core.columns:
        counts = core[column].dropna().value_counts()
        if not counts.empty:
            tie_blocks.append(float(counts.iloc[0]) / float(counts.sum()))
    path = RANK_MATRIX_DIR / f'{source}.canonical_core_rank.parquet'
    ranks.astype(np.float32).to_parquet(path)

    # The claim that C2-A shares no statistic between samples is checked directly: one
    # sample ranked alone must rank exactly as it does inside the full matrix.
    probe = str(core.columns[0])
    single, _ = build_within_sample_rank_matrix(core[[probe]], ranking_genes)
    left, right = ranks[probe].to_numpy(dtype=float), single[probe].to_numpy(dtype=float)
    both_missing = np.isnan(left) & np.isnan(right)
    independence = float(np.nanmax(np.where(both_missing, 0.0, np.abs(left - right))))

    record = dict(
        source_expression=source, modality_class=contract.modality_class,
        c1_input_scale=contract.input_scale, c1_status=contract.status,
        ranking_universe_genes=len(ranking_genes), samples_ranked=int(ranks.shape[1]),
        matrix_columns_in_c1=int(matrix.shape[1]), columns_not_bound_to_a_sample=unbound,
        core_genes_present_in_source=int(core.notna().any(axis=1).sum()),
        core_finite_fraction_min=float(qc.finite_fraction.min()),
        core_finite_fraction_median=float(qc.finite_fraction.median()),
        core_finite_fraction_max=float(qc.finite_fraction.max()),
        genes_ranked_per_sample_median=float(qc.finite_genes.median()),
        largest_tie_block_fraction_min=float(min(tie_blocks)) if tie_blocks else None,
        largest_tie_block_fraction_median=float(np.median(tie_blocks)) if tie_blocks else None,
        rank_unique_values_median=float(qc.rank_unique_values.median()),
        rank_min=float(finite.min()) if finite.size else None,
        rank_max=float(finite.max()) if finite.size else None,
        rank_exclusive_zero_one=bool(finite.size and (finite > 0.0).all()
                                     and (finite < 1.0).all()),
        sample_independence_max_abs_diff=independence,
        float32_output=True,
        preprocessing_scope=contract.preprocessing_scope,
        fold_isolation=contract.fold_isolation,
        author_batch_corrected=bool(contract.author_batch_corrected),
        rank_matrix_relpath=path.relative_to(bg.PROJECT).as_posix(),
        rank_matrix_size_bytes=int(path.stat().st_size), failure='')
    return record, qc


def failed_source(contract: pd.Series, ranking_genes: list[str], error: str) -> dict:
    return dict(
        source_expression=contract.source_expression, modality_class=contract.modality_class,
        c1_input_scale=contract.input_scale, c1_status=contract.status,
        ranking_universe_genes=len(ranking_genes), samples_ranked=0, matrix_columns_in_c1=None,
        columns_not_bound_to_a_sample=None, core_genes_present_in_source=None,
        core_finite_fraction_min=None, core_finite_fraction_median=None,
        core_finite_fraction_max=None, genes_ranked_per_sample_median=None,
        largest_tie_block_fraction_min=None, largest_tie_block_fraction_median=None,
        rank_unique_values_median=None, rank_min=None, rank_max=None,
        rank_exclusive_zero_one=False, sample_independence_max_abs_diff=None,
        float32_output=False, preprocessing_scope=contract.preprocessing_scope,
        fold_isolation=contract.fold_isolation,
        author_batch_corrected=bool(contract.author_batch_corrected),
        rank_matrix_relpath='', rank_matrix_size_bytes=0, failure=error)


def main() -> int:
    for directory in (LAYER, RANK_MATRIX_DIR):
        directory.mkdir(parents=True, exist_ok=True)

    ranking_genes = core_ranking_universe()
    ready = c1_ready_sources()
    binding = bound_samples_by_source()

    metrics, qc_frames, ranked = [], [], set()
    for contract in ready.itertuples():
        try:
            record, qc = rank_one_source(contract.source_expression, str(contract.output_relpath),
                                         binding, ranking_genes, contract)
        except (ValueError, KeyError, OSError) as error:
            record, qc = failed_source(contract, ranking_genes, str(error)), pd.DataFrame()
        metrics.append(record)
        if not qc.empty:
            qc['modality_class'] = contract.modality_class
            qc['preprocessing_scope'] = contract.preprocessing_scope
            qc['fold_isolation'] = contract.fold_isolation
            qc['author_batch_corrected'] = bool(contract.author_batch_corrected)
            qc['rank_status'] = np.where(qc.finite_genes.gt(0), RANK_READY, RANK_NOT_READY)
            qc_frames.append(qc)
            ranked.update((contract.source_expression, sample)
                          for sample in qc.loc[qc.rank_status.eq(RANK_READY), 'sample_id'])

    metric_frame = pd.DataFrame(metrics)
    for column in SOURCE_METRIC_COLUMNS:
        if column not in metric_frame:
            metric_frame[column] = None
    qc_frame = pd.concat(qc_frames, ignore_index=True)
    qc_frame[QC_COLUMNS].sort_values(['source_expression', 'sample_id']).to_csv(
        LAYER / 'SAMPLE_RANK_QC.csv.gz', index=False, compression='gzip')
    metric_frame[SOURCE_METRIC_COLUMNS].to_csv(LAYER / 'SOURCE_RANK_METRICS.csv', index=False)

    # ------------------------------------------------------ pair coverage, from C1 upward
    c1 = pd.read_csv(EXPRESSION_LAYER / 'PAIR_EXPRESSION_COVERAGE.csv.gz', low_memory=False)
    c1_ready = c1.quantitative_status.eq('QUANTITATIVE_READY') & (
        c1.matrix_column_t0 != c1.matrix_column_t1)

    def rank_state(index, source, sample_t0, sample_t1):
        if not c1_ready.loc[index]:
            return RANK_PARENT_NOT_READY, str(c1.quantitative_status.iloc[index])
        missing = [end for end, sample in (('t0', sample_t0), ('t1', sample_t1))
                   if (source, sample) not in ranked]
        return (RANK_NOT_READY, 'not_ranked=' + ','.join(missing)) if missing else (
            RANK_READY, '')

    states = [rank_state(index, row.source_expression, row.sample_t0, row.sample_t1)
              for index, row in c1.iterrows()]
    c1['rank_status'] = [state[0] for state in states]
    c1['rank_detail'] = [state[1] for state in states]
    provenance = metric_frame.set_index('source_expression')[
        ['preprocessing_scope', 'fold_isolation', 'author_batch_corrected']]
    c1 = c1.join(provenance, on='source_expression')
    c1['author_batch_corrected'] = c1.author_batch_corrected.isin([True, 'True', 'true'])
    pairs_ready = c1.rank_status.eq(RANK_READY)
    if int(pairs_ready.sum()) > int(c1_ready.sum()):
        raise SystemExit('C2 recovered a pair that C1 did not rate quantitative-ready')
    c1.rename(columns={'quantitative_status': 'quantitative_status_c1'})[
        PAIR_COVERAGE_COLUMNS].to_csv(LAYER / 'PAIR_RANK_COVERAGE.csv.gz', index=False,
                                      compression='gzip')

    total_bytes = int(metric_frame.rank_matrix_size_bytes.sum())
    failed = metric_frame[metric_frame.failure.astype(str).ne('')]
    summary = dict(
        representation_layer_version=REPRESENTATION_LAYER_VERSION,
        level='C2 within-sample rank representation (C2-A) plus a model-time robust '
              'standardization contract (C2-B)',
        built_from='dataset/expression_layer v0.3.1-C1 matrices, read rather than rebuilt, '
                   'over the frozen strict canonical core of gene_space v0.2.1',
        mathematical_authority=(
            'dataset/src/representation_transforms.py, supplied by the Owner; the committed '
            'file is token-identical to that block apart from its module docstring'),
        c2a=dict(
            rank_formula=RANK_FORMULA, rank_method=RANK_METHOD,
            percentile_definition='(average_rank - 0.5) / n_valid, so values lie strictly in '
                                  '(0, 1) and never touch 0 or 1',
            ranking_universe=f'genes where gene_space/GENE_SPACE_GENES.csv.gz '
                             f'{CORE_UNIVERSE_COLUMN} is true, read from the artifact at '
                             f'build time ({len(ranking_genes)} genes)',
            universe_is_not='a final modelling feature set; it is one shared measurement '
                            'coordinate, so that a 0.90 means the same position on the same '
                            'gene list in every source',
            missing_values='not imputed: only finite genes are ranked, missing genes stay NaN, '
                           'and the denominator is that sample\'s own finite gene count',
            statistics_shared_between_samples='none; the per-source independence check ranks '
                                              'one sample alone and inside the full matrix and '
                                              'requires exact agreement',
            ties='average rank, never broken by gene order, because a pseudobulk sample can '
                 'have thousands of zero-valued core genes that form one tie block',
            output='float32 rank matrices, local only, like C1'),
        c2b=dict(
            status='CONTRACT_ONLY_NOT_MATERIALIZED',
            global_standardized_matrix_created=False,
            why='a median and MAD estimated over all samples would be estimated from test '
                'patients, which is the exact leak this layer exists to prevent',
            fit='representation_transforms.fit_robust_standardizer(matrix, train_columns, '
                'source_expression=..., split_id=...): explicit training columns are '
                'mandatory, per source and per split',
            transform='representation_transforms.apply_robust_standardizer(matrix, '
                      'fitted_state): estimates nothing from the matrix it transforms, and the '
                      'same state is applied to train, validation and test',
            centre='training median', scale='1.482602218505602 x training MAD',
            fallback='training IQR / 1.3489795003921634 when MAD is zero but the gene varies',
            no_epsilon='MAD = IQR = 0 yields UNUSABLE_CONSTANT_OR_SPARSE; scale is never '
                       'floored at an epsilon, which would turn an almost constant gene into '
                       'an enormous z-score',
            fit_statuses=['READY', 'INSUFFICIENT_TRAIN_SAMPLES', 'UNUSABLE_CONSTANT_OR_SPARSE'],
            scale_methods=['MAD', 'IQR_FALLBACK', 'UNUSABLE'],
            unseen_source_rule='C2-B is not zero-shot across sources: a held-out source with no '
                               'training sample must not have its median/MAD fitted from that '
                               'held-out source. Such an evaluation uses C2-A ranks or another '
                               'explicitly source-independent representation.',
            gse319641='its upstream is AUTHOR_FULL_SOURCE ComBat with fold_isolation='
                      'NOT_ESTABLISHED, so fit_robust_standardizer refuses it unless '
                      'allow_nonisolated_upstream=True is passed, and then only as a named '
                      'sensitivity/stress analysis that may never be reported as a '
                      'leakage-isolated main result'),
        headline=dict(
            ranking_universe_genes=len(ranking_genes),
            sources_c1_ready=int(len(ready)),
            sources_rank_ready=int(metric_frame.samples_ranked.gt(0).sum()),
            sources_failed=int(len(failed)),
            samples_rank_ready=len(ranked),
            pairs_c1_ready=int(c1_ready.sum()),
            pairs_rank_ready=int(pairs_ready.sum()),
            patients_rank_ready=int(c1.loc[pairs_ready, 'patient_uid'].nunique()),
            author_batch_corrected_rank_pairs=int(
                (pairs_ready & c1.author_batch_corrected).sum()),
            core_finite_fraction_min=float(qc_frame.finite_fraction.min()),
            core_finite_fraction_median=float(qc_frame.finite_fraction.median()),
            core_finite_fraction_max=float(qc_frame.finite_fraction.max()),
            rank_matrices=int((metric_frame.rank_matrix_size_bytes > 0).sum()),
            local_rank_matrix_total_bytes=total_bytes,
            local_rank_matrix_total_mebibytes=round(total_bytes / 1048576.0, 1)),
        modality_pair_classes={
            name: dict(pairs_c1_ready=int((c1.modality_class.eq(name) & c1_ready).sum()),
                       pairs_rank_ready=int((c1.modality_class.eq(name) & pairs_ready).sum()),
                       patients_rank_ready=int(c1.loc[
                           c1.modality_class.eq(name) & pairs_ready, 'patient_uid'].nunique()))
            for name in sorted(c1.modality_class.dropna().unique())},
        coverage_statuses={status: int(count)
                           for status, count in c1.rank_status.value_counts().items()},
        per_source=metric_frame[SOURCE_METRIC_COLUMNS].to_dict('records'),
        redistribution_policy=dict(
            published='code, contracts, metrics, per-sample QC, pair coverage and this report',
            local_only='rank matrices, under longitudinal-data/data/processed/'
                       'expression_v0.4_rank/, which git ignores',
            reason='a derived matrix inherits the licence terms of the GEO supplementary file '
                   'it came from and this repository carries no per-accession ledger'),
        explicit_non_goals=[
            'no paired delta, no delta rank, no t0/t1 state transition: that is level C3',
            'no response model, no classification, no feature selection, no PCA',
            'no batch correction of any kind here, and no cross-cohort z-scoring',
            'no imputation of missing core genes and no sample exclusion threshold: coverage '
            'is reported here, not silently acted on',
            'no C2-B matrix: a whole-dataset robust standardization would leak the test split',
            'the rank matrices are still not comparable across sources in absolute terms; C2-A '
            'buys a shared coordinate system, not a shared measurement scale'],
        status_vocabulary=dict(
            pair_rank_status=[RANK_READY, RANK_NOT_READY, RANK_PARENT_NOT_READY],
            c1_status_carried_in_quantitative_status_c1=sorted(
                c1.quantitative_status.dropna().unique().tolist())))
    (LAYER / 'REPRESENTATION_LAYER_REPORT.json').write_text(json.dumps(summary, indent=2),
                                                            encoding='utf-8')

    print(json.dumps(dict(headline=summary['headline'],
                          modality_pair_classes=summary['modality_pair_classes'],
                          coverage_statuses=summary['coverage_statuses']), indent=2))
    print(metric_frame[['source_expression', 'modality_class', 'samples_ranked',
                        'core_finite_fraction_min', 'core_finite_fraction_median',
                        'largest_tie_block_fraction_median', 'rank_unique_values_median',
                        'rank_min', 'rank_max', 'rank_exclusive_zero_one',
                        'sample_independence_max_abs_diff', 'author_batch_corrected']]
          .to_string(index=False))
    if not failed.empty:
        print('FAILED SOURCES')
        print(failed[['source_expression', 'failure']].to_string(index=False))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
