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
import time
from pathlib import Path

import numpy as np
import pandas as pd

import build_gene_space as bg
from representation_transforms import (RANK_FORMULA, RANK_METHOD, RobustSplit,
                                       build_within_sample_rank_matrix, id_gap, id_gap_text)

DATASET_ROOT = bg.DATASET_ROOT
LAYER = DATASET_ROOT / 'representation_layer'
EXPRESSION_LAYER = DATASET_ROOT / 'expression_layer'
# Local only, like the C1 matrices.
RANK_MATRIX_DIR = bg.PROJECT / 'longitudinal-data' / 'data' / 'processed' / 'expression_v0.4_rank'

REPRESENTATION_LAYER_VERSION = 'v0.4.1-C2'
CORE_UNIVERSE_COLUMN = 'in_core_gene_space'
BUILD_COMPLETE = 'COMPLETE'
BUILD_INCOMPLETE = 'INCOMPLETE'

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
    'rank_exclusive_zero_one', 'float32_output', 'preprocessing_scope', 'fold_isolation',
    'author_batch_corrected', 'build_id', 'rank_matrix_relpath',
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


def rank_invariant_violation(ranks: np.ndarray, finite_expression_values: int) -> str:
    """'' when a rank matrix is a valid within-sample percentile set.

    C2-A's claims are that every measured gene gets a percentile and that the percentile
    grid is strictly inside (0, 1), so both properties stop the build rather than being
    reported and waved through.
    """
    finite = ranks[np.isfinite(ranks)]
    if not finite.size:
        return 'no finite rank value was produced'
    if finite.size != int(finite_expression_values):
        return (f'{int(finite_expression_values) - finite.size} finite gene value(s) received '
                'no rank: C2-A may neither impute what C1 lacked nor drop what C1 measured')
    if bool((finite <= 0.0).any()) or bool((finite >= 1.0).any()):
        return ('rank value touches 0 or 1: a midrank percentile is strictly inside (0, 1) '
                'by construction')
    return ''


def build_incompleteness(failed: pd.DataFrame) -> str:
    """Non-empty text means this build must not report success.

    Diagnostic artifacts are still written first, because a failed source is easier to fix
    from its metric row than from a traceback, but a partial build is not a build.
    """
    if failed.empty:
        return ''
    return (f'{len(failed)} QUANTITATIVE_READY C1 source(s) failed to rank: '
            + '; '.join(f'{row.source_expression}: {row.failure}'
                        for row in failed.head(5).itertuples()))


def membership_check(what: str, expected: list[str], actual: list[str]) -> str:
    """'' only when the two id sets are the same ids, not merely the same size.

    Used for the three memberships that define this layer: the genes a rank matrix holds,
    the samples it holds, and the intervals it rates ready.
    """
    text = id_gap_text(id_gap(expected, actual))
    return '' if not text else f'{what} does not match: {text}'


def rank_matrices_for_build(expected_build_id: str, layer_dir: Path | None = None) -> pd.DataFrame:
    """The only way a downstream level may obtain C2-A matrices: name the build.

    A build publishes only when it completed, so asking for an incomplete build, or for any
    build other than the one the layer currently reports, raises.  There is no fallback to
    older files: parquet left on disk by an earlier build is invisible to this function even
    though the path would still resolve, which is what makes a partial rebuild fail loudly
    instead of silently training on the last good matrices.
    """
    directory = LAYER if layer_dir is None else Path(layer_dir)
    report = json.loads((directory / 'REPRESENTATION_LAYER_REPORT.json').read_text('utf-8'))
    published = str(report['headline'].get('build_id') or '')
    status = str(report['headline'].get('build_status') or 'UNKNOWN')
    if status != BUILD_COMPLETE:
        raise SystemExit(f'{directory.name}: build {published or "unnamed"} is {status} and was '
                         'never published; refusing to serve its matrices')
    if expected_build_id != published:
        raise SystemExit(
            f'no published C2-A build {expected_build_id!r}; this layer publishes only '
            f'{published!r}, and it will not fall back to files left by an earlier build')
    metrics = pd.read_csv(directory / 'SOURCE_RANK_METRICS.csv')
    built = metrics[metrics.build_id.astype(str).eq(published)].copy()
    if built.empty:
        raise SystemExit(f'build {published!r} has a complete report but no metric rows')
    gap = membership_check('ranked sources', sorted(c1_ready_source_names()),
                           sorted(built.source_expression.astype(str)))
    if gap:
        raise SystemExit(f'build {published!r} is incomplete on disk: {gap}')
    return built[['source_expression', 'build_id', 'rank_matrix_relpath', 'samples_ranked',
                  'fold_isolation', 'author_batch_corrected']]


def c1_ready_source_names() -> set[str]:
    """Which sources C1 says are quantitative, read from C1 rather than restated."""
    metrics = pd.read_csv(EXPRESSION_LAYER / 'SOURCE_EXPRESSION_METRICS.csv')
    return set(metrics.loc[metrics.status.eq('QUANTITATIVE_READY'), 'source_expression'].astype(str))


def robust_split_for_source(source: str, split_id: str, train_sample_ids,
                            eval_sample_ids=()) -> RobustSplit:
    """Build a C2-B split whose fold isolation comes from the C1 manifest, not from a caller.

    The argument a caller must not get to choose is the one that decides whether a scaler
    may claim to be leakage-isolated: for GSE319641 the C1 layer recorded
    AUTHOR_FULL_SOURCE ComBat with fold_isolation NOT_ESTABLISHED, so the split built here
    carries that fact and fit_robust_standardizer refuses it.
    """
    metrics = pd.read_csv(EXPRESSION_LAYER / 'SOURCE_EXPRESSION_METRICS.csv')
    row = metrics[metrics.source_expression.eq(source)]
    if row.empty:
        raise ValueError(f'{source}: not a C1 source, so its fold isolation is unknown')
    if not row.status.eq('QUANTITATIVE_READY').all():
        raise ValueError(f'{source}: is not QUANTITATIVE_READY in C1')
    return RobustSplit(source_expression=source, split_id=split_id,
                       train_sample_ids=tuple(train_sample_ids),
                       eval_sample_ids=tuple(eval_sample_ids),
                       upstream_fold_isolation=str(row.fold_isolation.iloc[0]))


def rank_one_source(source: str, relpath: str, binding: pd.DataFrame,
                    ranking_genes: list[str], contract: pd.Series,
                    build_id: str) -> tuple[dict, pd.DataFrame]:
    """Rank one C1 source over the fixed core universe and write its float32 matrix locally."""
    matrix = pd.read_parquet(bg.PROJECT / relpath)
    bound = binding[binding.source_expression.eq(source)]
    labels = bound.set_index('matrix_column').sample_id
    columns = [column for column in matrix.columns if column in labels.index]
    unbound = int(matrix.shape[1] - len(columns))
    # C1 columns are native matrix headers; the sample identity lives in the binding
    # manifest.  This relabelling is 1:1 and in-memory only: the C1 file is never touched.
    named = matrix[columns].rename(columns=labels.to_dict())
    if named.columns.has_duplicates:
        raise ValueError(f'{source}: two native columns resolve to one sample id')

    core = named.reindex(pd.Index(ranking_genes))
    ranks, qc = build_within_sample_rank_matrix(core, ranking_genes)
    qc['source_expression'] = source

    values = ranks.to_numpy(dtype=float)
    finite = values[np.isfinite(values)]
    # The three identities this layer claims, checked by id rather than by count: the genes
    # are the frozen core, the samples are exactly what C1 bound for this source, and every
    # finite expression value carries a percentile strictly inside (0, 1).
    for problem in (
            membership_check(f'{source} gene rows', ranking_genes,
                             [str(gene) for gene in ranks.index]),
            membership_check(f'{source} sample columns',
                             sorted(bound.sample_id.astype(str)),
                             [str(sample) for sample in ranks.columns]),
            rank_invariant_violation(values,
                                     int(np.isfinite(core.to_numpy(dtype=float)).sum()))):
        if problem:
            raise ValueError(f'{source}: {problem}')

    tie_blocks = []
    for column in core.columns:
        counts = core[column].dropna().value_counts()
        if not counts.empty:
            tie_blocks.append(float(counts.iloc[0]) / float(counts.sum()))

    # Only a source that passed its own invariants reaches the disk.
    path = RANK_MATRIX_DIR / f'{source}.canonical_core_rank.parquet'
    ranks.astype(np.float32).to_parquet(path)

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
        float32_output=True,
        preprocessing_scope=contract.preprocessing_scope,
        fold_isolation=contract.fold_isolation,
        author_batch_corrected=bool(contract.author_batch_corrected),
        build_id=build_id,
        rank_matrix_relpath=path.relative_to(bg.PROJECT).as_posix(),
        rank_matrix_size_bytes=int(path.stat().st_size), failure='')
    return record, qc


def failed_source(contract: pd.Series, ranking_genes: list[str], error: str,
                  build_id: str) -> dict:
    return dict(
        source_expression=contract.source_expression, modality_class=contract.modality_class,
        c1_input_scale=contract.input_scale, c1_status=contract.status,
        ranking_universe_genes=len(ranking_genes), samples_ranked=0, matrix_columns_in_c1=None,
        columns_not_bound_to_a_sample=None, core_genes_present_in_source=None,
        core_finite_fraction_min=None, core_finite_fraction_median=None,
        core_finite_fraction_max=None, genes_ranked_per_sample_median=None,
        largest_tie_block_fraction_min=None, largest_tie_block_fraction_median=None,
        rank_unique_values_median=None, rank_min=None, rank_max=None,
        rank_exclusive_zero_one=False,
        float32_output=False, preprocessing_scope=contract.preprocessing_scope,
        fold_isolation=contract.fold_isolation,
        author_batch_corrected=bool(contract.author_batch_corrected),
        build_id=build_id,
        rank_matrix_relpath='', rank_matrix_size_bytes=0, failure=error)


def main() -> int:
    for directory in (LAYER, RANK_MATRIX_DIR):
        directory.mkdir(parents=True, exist_ok=True)

    ranking_genes = core_ranking_universe()
    ready = c1_ready_sources()
    binding = bound_samples_by_source()
    build_id = (f'{REPRESENTATION_LAYER_VERSION}+core{len(ranking_genes)}'
                f'+{time.strftime("%Y%m%dT%H%M%SZ", time.gmtime())}')

    metrics, qc_frames, ranked = [], [], set()
    for contract in ready.itertuples():
        try:
            record, qc = rank_one_source(contract.source_expression, str(contract.output_relpath),
                                         binding, ranking_genes, contract, build_id)
        except (ValueError, KeyError, OSError) as error:
            record, qc = failed_source(contract, ranking_genes, str(error), build_id), \
                pd.DataFrame()
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
    # A build in which every source failed still has to produce its diagnosis and exit
    # non-zero; concat([]) would crash before either happened.
    qc_frame = (pd.concat(qc_frames, ignore_index=True) if qc_frames
                else pd.DataFrame(columns=QC_COLUMNS))

    def coverage_statistic(statistic):
        values = pd.to_numeric(qc_frame.get('finite_fraction'), errors='coerce').dropna()
        return float(getattr(values, statistic)()) if len(values) else None

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
    recovered = sorted(set(c1.loc[pairs_ready, 'pair_uid'].astype(str))
                       - set(c1.loc[c1_ready, 'pair_uid'].astype(str)))
    if recovered:
        raise SystemExit(f'C2 rated {len(recovered)} interval(s) ready that C1 did not: '
                         + ','.join(recovered[:5]))
    c1.rename(columns={'quantitative_status': 'quantitative_status_c1'})[
        PAIR_COVERAGE_COLUMNS].to_csv(LAYER / 'PAIR_RANK_COVERAGE.csv.gz', index=False,
                                      compression='gzip')

    total_bytes = int(metric_frame.rank_matrix_size_bytes.sum())
    failed = metric_frame[metric_frame.failure.astype(str).ne('')]
    incomplete = build_incompleteness(failed)
    # A build publishes exactly once, when it completed; an incomplete build writes its
    # diagnosis and stays unpublished, so a consumer asking for it is refused.
    build_status = BUILD_INCOMPLETE if incomplete else BUILD_COMPLETE
    summary = dict(
        representation_layer_version=REPRESENTATION_LAYER_VERSION,
        level='C2 within-sample rank representation (C2-A) plus a model-time robust '
              'standardization contract (C2-B)',
        built_from='dataset/expression_layer v0.3.1-C1 matrices, read rather than rebuilt, '
                   'over the frozen strict canonical core of gene_space v0.2.1',
        mathematical_authority=(
            'dataset/src/representation_transforms.py, authored by the Owner: the C2-A rank '
            'and the C2-B estimator are his v0.4-C2 block, restructured on his instruction so '
            'that isolation and membership are enforced by what the code can see rather than '
            'by what a caller declares (see c2a.structural_guarantees and c2b.api_hardening)'),
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
            statistics_shared_between_samples='none, by structure: rank_one_sample receives '
                                              'one sample\'s vector and the matrix builder '
                                              'calls it once per column, so another sample is '
                                              'not visible in the ranking path',
            ties='average rank, never broken by gene order, because a pseudobulk sample can '
                 'have thousands of zero-valued core genes that form one tie block',
            structural_guarantees=[
                'sample isolation: the only ranking primitive takes a single pd.Series, and a '
                'test asserts rank_one_sample(column) == build_within_sample_rank_matrix(matrix)'
                '[column] for every column, with the spy recording that no call ever received '
                'more than one sample',
                'membership by identity: rank matrix gene rows must equal the frozen core and '
                'sample columns must equal the source\'s C1-bound sample ids, reported as '
                'missing/unexpected/duplicated ids rather than as counts',
                'value contract: a finite C1 gene value that receives no rank, or a rank that '
                'touches 0 or 1, fails the source before anything is written '
                '(rank_invariant_violation)'],
            output='float32 rank matrices, local only, like C1'),
        c2b=dict(
            status='CONTRACT_ONLY_NOT_MATERIALIZED',
            global_standardized_matrix_created=False,
            why='a median and MAD estimated over all samples would be estimated from test '
                'patients, which is the exact leak this layer exists to prevent',
            fit='representation_transforms.fit_robust_standardizer(matrix, split): the split '
                'is a RobustSplit carrying real train and eval membership, and the state that '
                'comes back records the samples and genes the fit actually touched',
            transform='representation_transforms.apply_robust_standardizer(matrix, state, '
                      'split): estimates nothing from the matrix it transforms, and the same '
                      'state is applied to train, validation and test',
            api_hardening=[
                'leakage is proved by membership, not by a declaration: '
                'verify_fit_membership(state, split) requires fit_sample_ids == '
                'split.train_sample_ids, and a split cannot be built with a sample on both '
                'sides, so a scaler fitted on more than the training membership cannot be '
                'applied at all - fitting on the full matrix fails on its own record',
                'the state carries its source and split and both must equal the split it is '
                'being applied to: a scaler from another cohort or another fold raises',
                'genes to transform must be exactly the genes that were fitted, and every '
                'sample column must belong to this split\'s train or eval side; duplicate gene '
                'rows and duplicate sample columns are rejected, since a repeated column '
                'would weight one sample twice',
                'upstream fold isolation travels on the split as data loaded from the C1 '
                'manifest by robust_split_for_source(); an empty or unknown value fails '
                'closed, so no caller inherits the safe answer',
                'the builder exits non-zero when any ready C1 source fails, and only a '
                'complete build is published to consumers'],
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
                      'NOT_ESTABLISHED, which robust_split_for_source copies off the C1 '
                      'manifest, so fit_robust_standardizer refuses it unless '
                      'allow_nonisolated_upstream=True is passed, and then only as a named '
                      'sensitivity/stress analysis that may never be reported as a '
                      'leakage-isolated main result'),
        headline=dict(
            build_id=build_id,
            build_status=build_status,
            published=build_status == BUILD_COMPLETE,
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
            core_finite_fraction_min=coverage_statistic('min'),
            core_finite_fraction_median=coverage_statistic('median'),
            core_finite_fraction_max=coverage_statistic('max'),
            rank_matrices=int((metric_frame.rank_matrix_size_bytes > 0).sum()),
            local_rank_matrix_total_bytes=total_bytes,
            local_rank_matrix_total_mebibytes=round(total_bytes / 1048576.0, 1)),
        matrix_consumption_rule=dict(
            loader='build_representation_layer.rank_matrices_for_build(expected_build_id), which '
                   'refuses a build that is not the published one and a build that did not '
                   'complete',
            publish_rule='a complete build is published; a partial build writes its diagnosis '
                         'and is never published, so asking for it raises',
            no_fallback='older matrices may remain on disk, and they are never served to a '
                        'consumer that named a different build',
            forbidden=['globbing expression_v0.4_rank/', 'taking the newest parquet in the '
                                                        'directory', 'reading a path that is not '
                                                                     'in SOURCE_RANK_METRICS.csv'],
            identity='build_id is published per source row and build_status per build in this '
                     'report; a source row says whether it ranked by leaving failure empty. '
                     'Identity is readable information, not a hash'),
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
                        'author_batch_corrected']]
          .to_string(index=False))
    if not failed.empty:
        print('FAILED SOURCES')
        print(failed[['source_expression', 'failure']].to_string(index=False))
    # Artifacts are written and diagnosable first; then a partial build refuses to look
    # like a complete one.
    incomplete = build_incompleteness(failed)
    if incomplete:
        raise SystemExit(f'C2 build incomplete: {incomplete}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
