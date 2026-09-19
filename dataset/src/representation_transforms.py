"""Level C2 representation mathematics, authored by the Owner.

The C2-A rank definition, the C2-B estimator and every rule below come from the Owner's
v0.4-C2 block unchanged.  The only later edit is the v0.4.1-C2 API hardening the Owner
ordered after reviewing that commit: upstream fold isolation must be declared rather than
defaulted, and a fitted state must name the source and split it belongs to.  Both are
fail-closed behaviour, not mathematics.

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


def midrank_percentile(
    values: pd.Series,
) -> pd.Series:
    """
    Convert one sample's gene-expression vector to within-sample
    percentile ranks.

    Only finite genes participate.

    Ties receive their average rank.

    Formula:
        percentile = (average_rank - 0.5) / n_valid

    Therefore values lie strictly inside (0, 1) when n_valid > 1.

    This function uses NO information from any other sample.
    """

    x = pd.to_numeric(values, errors="coerce").astype(float)

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
        Ranking is performed independently for every sample.
        No cohort/sample distribution is estimated.
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

        ranked = midrank_percentile(values)

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


def _validate_train_columns(
    matrix: pd.DataFrame,
    train_columns: Iterable[str],
) -> list[str]:

    requested = [str(c) for c in train_columns]

    if not requested:
        raise ValueError(
            "C2-B requires explicit non-empty training sample columns"
        )

    if len(set(requested)) != len(requested):
        raise ValueError("training sample columns contain duplicates")

    missing = [
        column
        for column in requested
        if column not in matrix.columns
    ]

    if missing:
        raise ValueError(
            "training samples are absent from source matrix: "
            + ", ".join(missing[:10])
        )

    return requested


def fit_robust_standardizer(
    canonical_matrix: pd.DataFrame,
    train_columns: Iterable[str],
    *,
    source_expression: str,
    split_id: str,
    upstream_fold_isolation: str | None = None,
    config: RobustFitConfig | None = None,
    allow_nonisolated_upstream: bool = False,
) -> pd.DataFrame:
    """
    Fit per-gene robust location/scale using TRAINING SAMPLES ONLY.

    Primary scale:
        1.4826 * MAD

    Fallback:
        IQR / 1.349

    If both are zero/non-estimable:
        gene is marked UNUSABLE_CONSTANT_OR_SPARSE

    No parameter is estimated from validation/test samples.

    upstream_fold_isolation has NO safe default on purpose: a caller that
    forgets it gets an error, not a favourable answer.
    """

    if config is None:
        config = RobustFitConfig()

    if upstream_fold_isolation is None:
        raise ValueError(
            f"{source_expression}: upstream_fold_isolation must be explicitly "
            "declared; the fold-safety of the representation this scaler is built "
            "from is not something a call may omit."
        )

    if (
        upstream_fold_isolation != "FOLD_INDEPENDENT"
        and not allow_nonisolated_upstream
    ):
        raise ValueError(
            f"{source_expression}: upstream representation has "
            f"fold_isolation={upstream_fold_isolation}; "
            "C2-B refuses to present this as leakage-isolated. "
            "Pass allow_nonisolated_upstream=True only for an explicitly "
            "flagged sensitivity/stress analysis."
        )

    train_columns = _validate_train_columns(
        canonical_matrix,
        train_columns,
    )

    train = (
        canonical_matrix[train_columns]
        .apply(pd.to_numeric, errors="coerce")
        .astype(float)
    )

    rows: list[dict] = []

    for gene, values in train.iterrows():

        arr = values.to_numpy(dtype=float)
        arr = arr[np.isfinite(arr)]

        n = len(arr)

        record = {
            "source_expression": source_expression,
            "split_id": split_id,
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

    state = pd.DataFrame(rows)

    if state.hgnc_symbol.duplicated().any():
        raise RuntimeError(
            "robust-standardizer state contains duplicate genes"
        )

    return state


def apply_robust_standardizer(
    canonical_matrix: pd.DataFrame,
    fitted_state: pd.DataFrame,
    *,
    expected_source_expression: str,
    expected_split_id: str,
) -> pd.DataFrame:
    """
    Apply a previously fitted C2-B transformer.

    This function estimates NOTHING from the matrix being transformed.

    The same fitted state can be applied to:
        train
        validation
        test

    Genes that were not estimable from the train fold remain NaN.

    The state carries the source and split it was fitted on, and the caller
    must name both: a scaler from another source, another split, or a state
    that mixes either is a different transformation, and applying it silently
    would be indistinguishable from a bug in a downstream pipeline.
    """

    required = {
        "hgnc_symbol",
        "median_train",
        "robust_scale",
        "fit_status",
        "source_expression",
        "split_id",
    }

    missing = required - set(fitted_state.columns)

    if missing:
        raise ValueError(
            f"fitted robust state missing columns: {sorted(missing)}"
        )

    for column, expected in (('source_expression', expected_source_expression),
                             ('split_id', expected_split_id)):
        observed = fitted_state[column].astype(str).unique().tolist()
        if len(observed) != 1:
            raise ValueError(
                f"fitted robust state mixes {len(observed)} values of {column}: "
                f"{sorted(observed)[:5]}; one state describes exactly one "
                f"{column}"
            )
        if observed[0] != str(expected):
            raise ValueError(
                f"fitted robust state was fitted on {column}={observed[0]}, not "
                f"{column}={expected}; refusing to apply a scaler across "
                f"{column} boundaries"
            )

    if fitted_state.hgnc_symbol.duplicated().any():
        raise ValueError("fitted robust state contains duplicate gene rows")

    if canonical_matrix.index.has_duplicates:
        raise ValueError("matrix to transform contains duplicate gene rows")

    if canonical_matrix.columns.has_duplicates:
        raise ValueError(
            "matrix to transform contains duplicate sample columns, which would "
            "weight one sample twice"
        )

    state = fitted_state.set_index("hgnc_symbol")

    genes = canonical_matrix.index.intersection(state.index)

    x = (
        canonical_matrix.loc[genes]
        .apply(pd.to_numeric, errors="coerce")
        .astype(float)
    )

    location = state.loc[genes, "median_train"]
    scale = state.loc[genes, "robust_scale"]
    status = state.loc[genes, "fit_status"]

    usable = status.eq("READY") & np.isfinite(scale) & scale.gt(0)

    output = pd.DataFrame(
        np.nan,
        index=canonical_matrix.index,
        columns=canonical_matrix.columns,
        dtype=np.float32,
    )

    usable_genes = genes[usable.to_numpy()]

    if len(usable_genes):
        standardized = (
            x.loc[usable_genes]
            .sub(location.loc[usable_genes], axis=0)
            .div(scale.loc[usable_genes], axis=0)
        )

        output.loc[usable_genes] = standardized.astype(np.float32)

    output.index.name = canonical_matrix.index.name or "hgnc_symbol"

    return output


def robust_fit_summary(
    fitted_state: pd.DataFrame,
) -> dict:
    """
    Compact fold-level QC for C2-B.
    """

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
