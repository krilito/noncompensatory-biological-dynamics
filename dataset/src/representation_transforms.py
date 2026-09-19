"""Level C2 representation mathematics, authored by the Owner.

The C2-A rank definition and the C2-B estimator are his v0.4-C2 block, and their
mathematics have not been edited.  What v0.4.1-C2 moved is where the guarantees live:
C2-A ranks come from rank_one_sample, which receives one sample's vector and therefore
cannot read another sample, and C2-B carries the samples and genes a fit actually
touched inside the state object.  Isolation and identity are shown by comparing
membership, not by a tolerance, a threshold, a hash or a caller's declaration.

C2-A (within-sample canonical-core percentile) is a fixed dataset representation; C2-B
(train-fold within-source robust standardization) is a model-time fit/transform contract and
must never be materialized as a whole-dataset standardized matrix.  Do not replace these
rules with sklearn StandardScaler, scipy zscore, global zscore, quantile normalization or
ComBat.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np
import pandas as pd


# ============================================================
# C2-A
# WITHIN-SAMPLE CANONICAL-GENE PERCENTILE
# ============================================================

RANK_METHOD = "average"
RANK_FORMULA = "MIDRANK_PERCENTILE"


def rank_one_sample(
    gene_values: pd.Series,
) -> pd.Series:
    """
    Convert ONE sample's gene-expression vector to within-sample
    percentile ranks.

    This is the whole of C2-A's computation: its argument is a single
    sample's vector, so no other sample's value is reachable from here, and
    the sample-independence property is a fact about this signature rather
    than a number someone hopes holds.

    Only finite genes participate.

    Ties receive their average rank.

    Formula:
        percentile = (average_rank - 0.5) / n_valid

    Therefore values lie strictly inside (0, 1) when n_valid > 1.
    """

    x = pd.to_numeric(gene_values, errors="coerce").astype(float)

    finite = np.isfinite(x.to_numpy())
    out = pd.Series(np.nan, index=x.index, dtype=float)

    n_valid = int(finite.sum())

    if n_valid == 0:
        return out

    if n_valid == 1:
        out.loc[finite] = 0.5
        return out

    valid = x.loc[finite]

    ranks = valid.rank(
        method=RANK_METHOD,
        ascending=True,
    )

    out.loc[valid.index] = (
        ranks.to_numpy(dtype=float) - 0.5
    ) / float(n_valid)

    return out


def build_within_sample_rank_matrix(
    canonical_matrix: pd.DataFrame,
    ranking_genes: Iterable[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Build C2-A representation.

    Input:
        rows    = HGNC canonical genes
        columns = samples
        values  = C1 quantitative expression

    ranking_genes:
        fixed dataset-level canonical ranking universe.

        For PAIR C2-A this should be the frozen strict canonical core
        from gene-space v0.2.1 (currently 9,338 genes).

    Returns:
        rank_matrix:
            rows = ranking genes
            columns = samples
            values = within-sample percentile ranks

        sample_qc:
            one row per sample describing finite core-gene coverage.

    IMPORTANT:
        This function does no ranking of its own.  It calls rank_one_sample
        on one column at a time and assembles the results, so a sample's
        ranks cannot depend on any other sample except through this loop.
    """

    genes = pd.Index([str(g) for g in ranking_genes])

    if genes.has_duplicates:
        raise ValueError("ranking gene universe contains duplicates")

    if canonical_matrix.index.has_duplicates:
        raise ValueError(
            "canonical expression matrix contains duplicate gene rows"
        )

    x = canonical_matrix.reindex(genes)

    ranks = pd.DataFrame(
        index=genes,
        columns=x.columns,
        dtype=np.float32,
    )

    qc_rows: list[dict] = []

    total_genes = len(genes)

    for sample in x.columns:
        values = x[sample]

        finite_count = int(
            np.isfinite(
                pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
            ).sum()
        )

        finite_fraction = (
            finite_count / total_genes
            if total_genes
            else 0.0
        )

        ranked = rank_one_sample(values)

        ranks[sample] = ranked.astype(np.float32)

        finite_ranks = ranked[np.isfinite(ranked)]

        qc_rows.append(
            {
                "sample_id": str(sample),
                "ranking_universe_genes": total_genes,
                "finite_genes": finite_count,
                "finite_fraction": finite_fraction,
                "rank_min": (
                    float(finite_ranks.min())
                    if len(finite_ranks)
                    else np.nan
                ),
                "rank_max": (
                    float(finite_ranks.max())
                    if len(finite_ranks)
                    else np.nan
                ),
                "rank_unique_values": (
                    int(finite_ranks.nunique())
                    if len(finite_ranks)
                    else 0
                ),
            }
        )

    ranks.index.name = "hgnc_symbol"

    qc = pd.DataFrame(qc_rows)

    return ranks, qc


# ============================================================
# C2-B
# TRAIN-FOLD ROBUST STANDARDIZATION
# ============================================================

MAD_TO_SIGMA = 1.482602218505602
IQR_TO_SIGMA = 1.3489795003921634


@dataclass(frozen=True)
class RobustFitConfig:
    """
    C2-B is a MODEL-TIME transform.

    Parameters must be estimated only from training samples of one source.
    """

    min_train_samples_per_gene: int = 4

    # We do NOT add arbitrary epsilon to a zero scale.
    # If MAD and IQR are both zero, the gene is unusable in that fold.
    allow_iqr_fallback: bool = True


def id_gap(expected: Iterable[str], actual: Iterable[str]) -> dict[str, list[str]]:
    """
    Which ids are missing, unexpected or duplicated, by identity.

    A count comparison passes for a set that lost one gene and gained another, which is
    precisely the corruption a membership check exists to catch.
    """
    wanted = {str(value) for value in expected}
    got = [str(value) for value in actual]
    counts: dict[str, int] = {}
    for value in got:
        counts[value] = counts.get(value, 0) + 1
    return {
        'missing': sorted(wanted - set(counts)),
        'unexpected': sorted(set(counts) - wanted),
        'duplicated': sorted(value for value, count in counts.items() if count > 1),
    }


def id_gap_text(gap: dict[str, list[str]]) -> str:
    return '' if not any(gap.values()) else '; '.join(
        f'{kind}={len(ids)}[{",".join(ids[:5])}]' for kind, ids in sorted(gap.items()) if ids)


@dataclass(frozen=True)
class RobustSplit:
    """
    A real split definition: the membership of one source, not a claim about it.

    train_sample_ids and eval_sample_ids are disjoint by construction, so a scaler
    fitted on the train ids is proven not to have seen an evaluation sample by
    comparing two lists, which is a fact about data rather than about an argument.

    upstream_fold_isolation is carried here from the C1 manifest as data.  It is
    empty unless someone loaded it, and an empty value fails closed at fit time.
    """

    source_expression: str
    split_id: str
    train_sample_ids: tuple[str, ...]
    eval_sample_ids: tuple[str, ...] = ()
    upstream_fold_isolation: str = ''

    def __post_init__(self) -> None:
        object.__setattr__(self, 'source_expression', str(self.source_expression))
        object.__setattr__(self, 'split_id', str(self.split_id))
        object.__setattr__(self, 'train_sample_ids',
                           tuple(str(value) for value in self.train_sample_ids))
        object.__setattr__(self, 'eval_sample_ids',
                           tuple(str(value) for value in self.eval_sample_ids))
        object.__setattr__(self, 'upstream_fold_isolation',
                           str(self.upstream_fold_isolation or ''))
        if not self.source_expression or not self.split_id:
            raise ValueError('a split must name both its source and its split id')
        if not self.train_sample_ids:
            raise ValueError(
                f'{self.source_expression}: a split with no training sample cannot fit anything')
        for label, ids in (('train', self.train_sample_ids), ('eval', self.eval_sample_ids)):
            if len(set(ids)) != len(ids):
                raise ValueError(f'{self.source_expression}/{self.split_id}: {label} sample ids '
                                 'contain duplicates')
        overlap = set(self.train_sample_ids) & set(self.eval_sample_ids)
        if overlap:
            raise ValueError(
                f'{self.source_expression}/{self.split_id}: {len(overlap)} sample(s) are both '
                f'training and evaluation samples [{",".join(sorted(overlap)[:5])}]')

    @property
    def known_sample_ids(self) -> set[str]:
        return set(self.train_sample_ids) | set(self.eval_sample_ids)


@dataclass(frozen=True)
class RobustState:
    """
    A fitted C2-B transformer plus the membership that produced it.

    Identity and sample/gene sets travel with the numbers because a DataFrame of
    per-gene statistics cannot show what it was fitted on: the same 9,338 rows
    could have come from any fold of any cohort.
    """

    frame: pd.DataFrame
    source_expression: str
    split_id: str
    fit_sample_ids: tuple[str, ...]
    fit_gene_ids: tuple[str, ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, 'fit_sample_ids',
                           tuple(str(value) for value in self.fit_sample_ids))
        object.__setattr__(self, 'fit_gene_ids',
                           tuple(str(value) for value in self.fit_gene_ids))


def verify_fit_membership(
    state: RobustState,
    split: RobustSplit,
) -> str:
    """
    '' only when the state was fitted on exactly what this split trains on.

    This is the leakage proof: fit_sample_ids == train_sample_ids, and train is
    disjoint from eval by the split's own construction, so the fitted numbers
    cannot have come from an evaluation sample.
    """
    if state.source_expression != split.source_expression:
        return (f'scaler belongs to source {state.source_expression}, this split trains '
                f'{split.source_expression}')
    if state.split_id != split.split_id:
        return (f'scaler belongs to split {state.split_id}, this split is {split.split_id}')
    gap = id_gap(split.train_sample_ids, state.fit_sample_ids)
    text = id_gap_text(gap)
    if text:
        return ('scaler was not fitted on exactly this split\'s training samples: ' + text)
    if set(state.fit_sample_ids) & set(split.eval_sample_ids):
        return 'scaler was fitted on samples this split evaluates'
    return ''


def fit_robust_standardizer(
    canonical_matrix: pd.DataFrame,
    split: RobustSplit,
    *,
    config: RobustFitConfig | None = None,
    allow_nonisolated_upstream: bool = False,
) -> RobustState:
    """
    Fit per-gene robust location/scale using TRAINING SAMPLES ONLY.

    Primary scale:
        1.4826 * MAD

    Fallback:
        IQR / 1.349

    If both are zero/non-estimable:
        gene is marked UNUSABLE_CONSTANT_OR_SPARSE

    The returned state records which samples and genes it actually touched,
    and fit verifies that record against the split.  Fitting on more than the
    training membership is therefore not something a caller may declare away:
    it is a contradiction the function detects in its own output.
    """

    if config is None:
        config = RobustFitConfig()

    if (
        split.upstream_fold_isolation != "FOLD_INDEPENDENT"
        and not allow_nonisolated_upstream
    ):
        raise ValueError(
            f"{split.source_expression}: upstream representation has "
            f"fold_isolation={split.upstream_fold_isolation or 'UNDECLARED'}; "
            "C2-B refuses to present this as leakage-isolated. "
            "Pass allow_nonisolated_upstream=True only for an explicitly "
            "flagged sensitivity/stress analysis."
        )

    if canonical_matrix.index.has_duplicates:
        raise ValueError(
            f'{split.source_expression}: matrix to fit contains duplicate gene rows')

    if canonical_matrix.columns.has_duplicates:
        raise ValueError(
            f'{split.source_expression}: matrix to fit contains duplicate sample columns, '
            'which would weight one sample twice')

    gap = id_gap(split.train_sample_ids, canonical_matrix.columns)
    if any(gap.values()):
        # A fit may not be handed the whole source and told to ignore the rest.  Estimating
        # only from training samples is a property of what this function was allowed to
        # read, so anything outside the training membership - including an evaluation
        # sample of this very split - stops the fit rather than being sliced away.
        raise ValueError(
            f'{split.source_expression}/{split.split_id}: a robust fit may only be handed its '
            'own training samples, but the matrix given is not exactly that membership: '
            + (id_gap_text(gap) or 'no ids to compare'))

    train = (
        canonical_matrix[list(split.train_sample_ids)]
        .apply(pd.to_numeric, errors="coerce")
        .astype(float)
    )

    rows: list[dict] = []

    for gene, values in train.iterrows():

        arr = values.to_numpy(dtype=float)
        arr = arr[np.isfinite(arr)]

        n = len(arr)

        record = {
            "source_expression": split.source_expression,
            "split_id": split.split_id,
            "hgnc_symbol": str(gene),
            "n_train_finite": n,
            "median_train": np.nan,
            "mad_train": np.nan,
            "iqr_train": np.nan,
            "robust_scale": np.nan,
            "scale_method": "UNUSABLE",
            "fit_status": "",
        }

        if n < config.min_train_samples_per_gene:
            record["fit_status"] = "INSUFFICIENT_TRAIN_SAMPLES"
            rows.append(record)
            continue

        median = float(np.median(arr))
        abs_dev = np.abs(arr - median)
        mad = float(np.median(abs_dev))

        q25, q75 = np.percentile(
            arr,
            [25.0, 75.0],
            method="linear",
        )

        iqr = float(q75 - q25)

        record["median_train"] = median
        record["mad_train"] = mad
        record["iqr_train"] = iqr

        mad_scale = MAD_TO_SIGMA * mad

        if np.isfinite(mad_scale) and mad_scale > 0.0:
            record["robust_scale"] = mad_scale
            record["scale_method"] = "MAD"
            record["fit_status"] = "READY"

        elif config.allow_iqr_fallback:
            iqr_scale = iqr / IQR_TO_SIGMA

            if np.isfinite(iqr_scale) and iqr_scale > 0.0:
                record["robust_scale"] = iqr_scale
                record["scale_method"] = "IQR_FALLBACK"
                record["fit_status"] = "READY"
            else:
                record["fit_status"] = (
                    "UNUSABLE_CONSTANT_OR_SPARSE"
                )

        else:
            record["fit_status"] = (
                "UNUSABLE_CONSTANT_OR_SPARSE"
            )

        rows.append(record)

    frame = pd.DataFrame(rows)

    if frame.hgnc_symbol.duplicated().any():
        raise RuntimeError(
            "robust-standardizer state contains duplicate genes"
        )

    # The record is read back from what the fit actually touched, not restated from the
    # request, so an implementation that fitted more than the training membership is
    # caught by its own output.
    state = RobustState(
        frame=frame,
        source_expression=split.source_expression,
        split_id=split.split_id,
        fit_sample_ids=tuple(str(column) for column in train.columns),
        fit_gene_ids=tuple(str(gene) for gene in frame.hgnc_symbol),
    )

    problem = verify_fit_membership(state, split)

    if problem:
        raise ValueError(
            f'{split.source_expression}/{split.split_id}: {problem}'
        )

    return state


def apply_robust_standardizer(
    canonical_matrix: pd.DataFrame,
    fitted_state: RobustState,
    split: RobustSplit,
) -> pd.DataFrame:
    """
    Apply a previously fitted C2-B transformer.

    This function estimates NOTHING from the matrix being transformed.

    The same fitted state can be applied to:
        train
        validation
        test

    Genes that were not estimable from the train fold remain NaN.

    What makes this safe is membership, not a declaration: the state must have
    been fitted on exactly this split's training samples, the split's train and
    eval sides are disjoint by construction, and every sample or gene in the
    matrix must be one this split already accounts for.
    """

    problem = verify_fit_membership(fitted_state, split)

    if problem:
        raise ValueError(
            f'{split.source_expression}/{split.split_id}: {problem}'
        )

    if fitted_state.frame.hgnc_symbol.duplicated().any():
        raise ValueError("fitted robust state contains duplicate gene rows")

    if canonical_matrix.index.has_duplicates:
        raise ValueError("matrix to transform contains duplicate gene rows")

    if canonical_matrix.columns.has_duplicates:
        raise ValueError(
            "matrix to transform contains duplicate sample columns, which would "
            "weight one sample twice"
        )

    names = [str(value) for value in canonical_matrix.index]

    gene_text = id_gap_text(id_gap(fitted_state.fit_gene_ids, names))

    if gene_text:
        raise ValueError(
            f'{split.source_expression}/{split.split_id}: the genes to transform are not the '
            'genes that were fitted: ' + gene_text
        )

    foreign = sorted(set(str(column) for column in canonical_matrix.columns)
                     - split.known_sample_ids)

    if foreign:
        raise ValueError(
            f'{split.source_expression}/{split.split_id}: {len(foreign)} sample(s) belong to '
            f'neither this split\'s training nor its evaluation side: '
            f'{",".join(foreign[:5])}'
        )

    state = fitted_state.frame.set_index("hgnc_symbol").loc[names]

    x = (
        canonical_matrix
        .apply(pd.to_numeric, errors="coerce")
        .astype(float)
    )

    location = state["median_train"].to_numpy(dtype=float)[:, None]
    scale = state["robust_scale"].to_numpy(dtype=float)[:, None]
    usable = (
        state["fit_status"].eq("READY").to_numpy()
        & np.isfinite(scale[:, 0])
        & (scale[:, 0] > 0.0)
    )[:, None]

    with np.errstate(divide="ignore", invalid="ignore"):
        standardized = np.where(usable, (x.to_numpy(dtype=float) - location) / scale, np.nan)

    output = pd.DataFrame(
        standardized.astype(np.float32),
        index=canonical_matrix.index,
        columns=canonical_matrix.columns,
        dtype=np.float32,
    )

    output.index.name = canonical_matrix.index.name or "hgnc_symbol"

    return output


def robust_fit_summary(
    fitted_state,
) -> dict:
    """
    Compact fold-level QC for C2-B.
    """

    fitted_state = fitted_state.frame if isinstance(fitted_state, RobustState) else fitted_state

    total = len(fitted_state)

    ready = fitted_state.fit_status.eq("READY")

    return {
        "genes_total": total,
        "genes_ready": int(ready.sum()),
        "genes_ready_fraction": (
            float(ready.mean()) if total else 0.0
        ),
        "genes_using_mad": int(
            fitted_state.scale_method.eq("MAD").sum()
        ),
        "genes_using_iqr_fallback": int(
            fitted_state.scale_method.eq("IQR_FALLBACK").sum()
        ),
        "genes_insufficient_train_samples": int(
            fitted_state.fit_status.eq(
                "INSUFFICIENT_TRAIN_SAMPLES"
            ).sum()
        ),
        "genes_constant_or_sparse": int(
            fitted_state.fit_status.eq(
                "UNUSABLE_CONSTANT_OR_SPARSE"
            ).sum()
        ),
    }
