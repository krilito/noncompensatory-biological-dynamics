"""Level C1 mathematical authority for the PAIR expression layer.

The transformation rules below are the contract supplied for Expression Layer
v0.3-C1.  They are deliberately NOT a generic normalization pipeline: each
declared source scale has exactly one legal route to a canonical-gene matrix, and
a source whose semantics cannot be established fails closed.

Two scale members were added to the supplied contract because five locally present
sources do not fit any supplied member.  Both additions are additive: they reuse
the existing ``log2(x + 1)`` -> median-per-gene route, they are marked in the
metrics, and neither performs any cross-cohort or library-size renormalization.
See ``dataset/expression_layer/README.md`` section
"Declared-scale extensions" for the adjudication that motivated them.
"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Literal

import numpy as np
import pandas as pd


class ExpressionScale(str, Enum):
    # Additive count scales
    RAW_COUNTS = "RAW_COUNTS"

    # Already library-normalized but still linear
    TPM = "TPM"
    FPKM = "FPKM"

    # Already transformed quantitative scales
    LOG2_CPM = "LOG2_CPM"
    LOG2_TPM = "LOG2_TPM"
    LOG2_FPKM = "LOG2_FPKM"
    LOG_EXPRESSION = "LOG_EXPRESSION"

    # Microarray / array source already normalized and log-like
    ARRAY_NORMALIZED_LOG = "ARRAY_NORMALIZED_LOG"

    # EXTENSION 1: array source normalized by the platform's own procedure but
    # still on a LINEAR intensity scale (MAS5, Illumina BASE quantile).  It is not
    # log-like, so it must not take the identity route, and it must not be
    # renormalized by library size, because array intensity is not counts.
    ARRAY_NORMALIZED_LINEAR = "ARRAY_NORMALIZED_LINEAR"

    # EXTENSION 2: linear-scale abundance matrix that the AUTHOR already
    # batch-corrected.  Provenance must state the exact chain, because the route
    # below recovers the pre-back-transform log-scale state the author analysed.
    # Never pool such a source as if it were uncorrected.
    LINEAR_ABUNDANCE_COMBAT_ADJUSTED = "LINEAR_ABUNDANCE_COMBAT_ADJUSTED"

    # Do not guess unknown scales
    UNKNOWN = "UNKNOWN"


@dataclass(frozen=True)
class SourceExpressionContract:
    source_expression: str
    modality: str
    scale: ExpressionScale

    # Whether native features are probes rather than gene-level measurements.
    probe_based: bool = False

    # Must be supported by source metadata / paper / file format.
    normalization_provenance: str = ""

    # Aggregation policy is explicit, not inferred at runtime.
    gene_aggregation: Literal[
        "SUM_BEFORE_TRANSFORM",
        "MEDIAN_AFTER_TRANSFORM",
    ] = "MEDIAN_AFTER_TRANSFORM"


COUNT_SCALES = {
    ExpressionScale.RAW_COUNTS,
}

LINEAR_ABUNDANCE_SCALES = {
    ExpressionScale.TPM,
    ExpressionScale.FPKM,
}

ALREADY_TRANSFORMED_SCALES = {
    ExpressionScale.LOG2_CPM,
    ExpressionScale.LOG2_TPM,
    ExpressionScale.LOG2_FPKM,
    ExpressionScale.LOG_EXPRESSION,
    ExpressionScale.ARRAY_NORMALIZED_LOG,
}

# EXTENSION sets.  Both take log2(x + 1) and then the median per canonical gene.
LINEAR_INTENSITY_SCALES = {
    ExpressionScale.ARRAY_NORMALIZED_LINEAR,
}

AUTHOR_BATCH_CORRECTED_SCALES = {
    ExpressionScale.LINEAR_ABUNDANCE_COMBAT_ADJUSTED,
}

LOG2_OFFSET_SCALES = LINEAR_INTENSITY_SCALES | AUTHOR_BATCH_CORRECTED_SCALES


def _numeric_frame(values: pd.DataFrame) -> pd.DataFrame:
    out = values.apply(pd.to_numeric, errors="coerce")

    if out.isna().all(axis=None):
        raise ValueError("expression matrix contains no numeric values")

    return out


def validate_scale(
    values: pd.DataFrame,
    contract: SourceExpressionContract,
) -> dict:
    """
    Validate declared source semantics.

    IMPORTANT:
    This function does not auto-detect and replace the declared scale.
    It only checks whether the declared contract is plausible.
    """

    x = _numeric_frame(values)
    arr = x.to_numpy(dtype=float)

    finite = arr[np.isfinite(arr)]

    if finite.size == 0:
        raise ValueError(f"{contract.source_expression}: no finite values")

    negative_fraction = float(np.mean(finite < 0))
    zero_fraction = float(np.mean(finite == 0))

    result = {
        "min": float(np.min(finite)),
        "max": float(np.max(finite)),
        "negative_fraction": negative_fraction,
        "zero_fraction": zero_fraction,
    }

    if contract.scale in COUNT_SCALES:
        if negative_fraction > 0:
            raise ValueError(
                f"{contract.source_expression}: RAW_COUNTS contains negative values"
            )

        integer_like = np.isclose(finite, np.round(finite), atol=1e-8)
        integer_fraction = float(np.mean(integer_like))

        result["integer_like_fraction"] = integer_fraction

        # Fail closed. Do not silently treat arbitrary continuous values as counts.
        if integer_fraction < 0.99:
            raise ValueError(
                f"{contract.source_expression}: declared RAW_COUNTS but "
                f"only {integer_fraction:.3%} are integer-like"
            )

    elif contract.scale in LINEAR_ABUNDANCE_SCALES:
        if negative_fraction > 0:
            raise ValueError(
                f"{contract.source_expression}: {contract.scale.value} "
                "contains negative values"
            )

    elif contract.scale in ALREADY_TRANSFORMED_SCALES:
        # Negative values can be legitimate after centering/log preprocessing.
        pass

    elif contract.scale in LINEAR_INTENSITY_SCALES:
        # EXTENSION 1 validation: a linear intensity may be zero but never
        # negative, and it is not a count, so no integer contract applies.
        if negative_fraction > 0:
            raise ValueError(
                f"{contract.source_expression}: {contract.scale.value} "
                "contains negative values"
            )

    elif contract.scale in AUTHOR_BATCH_CORRECTED_SCALES:
        # EXTENSION 2 validation: the declared route is log2(x + 1), which is
        # undefined below -1.  A matrix whose documented chain is
        # 2 ** y - 1 must therefore respect that floor.
        if float(np.min(finite)) < -1.0:
            raise ValueError(
                f"{contract.source_expression}: declared "
                f"{contract.scale.value} but values fall below the -1 floor of "
                "the documented 2 ** y - 1 back-transform"
            )
        result["author_batch_corrected"] = True

    elif contract.scale == ExpressionScale.UNKNOWN:
        raise ValueError(
            f"{contract.source_expression}: expression scale is UNKNOWN; "
            "adjudicate source semantics before transformation"
        )

    else:
        raise ValueError(f"unsupported expression scale: {contract.scale}")

    return result


def log2_cpm(counts: pd.DataFrame) -> pd.DataFrame:
    """
    Deterministic per-sample CPM followed by log2(CPM + 1).

    Genes/features are rows.
    Samples are columns.
    """

    x = _numeric_frame(counts).astype(float)

    library_size = x.sum(axis=0)

    if (library_size <= 0).any():
        bad = list(library_size.index[library_size <= 0])
        raise ValueError(
            f"zero/non-positive library size for samples: {bad[:10]}"
        )

    cpm = x.divide(library_size, axis=1) * 1_000_000.0

    return np.log2(cpm + 1.0)


def transform_continuous(
    values: pd.DataFrame,
    scale: ExpressionScale,
) -> pd.DataFrame:
    """
    Transform non-count expression without inventing cross-platform comparability.
    """

    x = _numeric_frame(values).astype(float)

    if scale == ExpressionScale.TPM:
        if (x < 0).any(axis=None):
            raise ValueError("TPM contains negative values")
        return np.log2(x + 1.0)

    if scale == ExpressionScale.FPKM:
        if (x < 0).any(axis=None):
            raise ValueError("FPKM contains negative values")
        return np.log2(x + 1.0)

    if scale in ALREADY_TRANSFORMED_SCALES:
        # Never log a value twice.
        return x.copy()

    if scale in LOG2_OFFSET_SCALES:
        # EXTENSION: a linear quantity whose only legal within-source log
        # representation is log2(x + 1).  No library-size renormalization.
        return np.log2(x + 1.0)

    raise ValueError(
        f"transform_continuous called for unsupported scale: {scale}"
    )


def mapped_feature_table(
    native_values: pd.DataFrame,
    feature_map: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.Series]:
    """
    Restrict to unambiguous canonical-gene mappings.

    native_values:
        rows = original_feature_id
        columns = samples

    Returns:
        filtered values
        HGNC-symbol series aligned to rows
    """

    required = {
        "original_feature_id",
        "hgnc_symbol",
        "mapping_status",
    }

    missing = required - set(feature_map.columns)
    if missing:
        raise ValueError(f"feature map missing columns: {sorted(missing)}")

    fm = feature_map.copy()

    accepted = {
        "EXACT",
        "ALIAS",
        "PREVIOUS_SYMBOL",
        "PLATFORM_ANNOTATION",
    }

    fm = fm[
        fm["mapping_status"].isin(accepted)
        & fm["hgnc_symbol"].notna()
    ].copy()

    # One canonical assignment per native feature.
    fm = fm.drop_duplicates("original_feature_id", keep=False)

    common = native_values.index.intersection(fm["original_feature_id"])

    x = native_values.loc[common].copy()

    gene_lookup = (
        fm.set_index("original_feature_id")
        .loc[common, "hgnc_symbol"]
        .astype(str)
    )

    return x, gene_lookup


def aggregate_counts_to_genes(
    native_counts: pd.DataFrame,
    genes: pd.Series,
) -> pd.DataFrame:
    """
    Count-like measurements are additive.

    Multiple native rows resolving to the same canonical gene are summed
    BEFORE library-size normalization.
    """

    x = native_counts.copy()
    x["__gene__"] = genes.to_numpy()

    grouped = x.groupby("__gene__", sort=True).sum(numeric_only=True)

    grouped.index.name = "hgnc_symbol"

    return grouped


def aggregate_continuous_to_genes(
    transformed_values: pd.DataFrame,
    genes: pd.Series,
) -> pd.DataFrame:
    """
    Continuous / normalized / probe-level measurements are not additive.

    Multiple features mapping to one gene are collapsed using the median.
    This applies especially to microarray probes.
    """

    x = transformed_values.copy()
    x["__gene__"] = genes.to_numpy()

    grouped = x.groupby("__gene__", sort=True).median(numeric_only=True)

    grouped.index.name = "hgnc_symbol"

    return grouped


def build_canonical_quantitative_matrix(
    native_values: pd.DataFrame,
    feature_map: pd.DataFrame,
    contract: SourceExpressionContract,
) -> tuple[pd.DataFrame, dict]:
    """
    Main deterministic transformation.

    Output:
        canonical genes x samples

    IMPORTANT:
    This produces a source-valid quantitative representation.
    It does NOT make different cohorts/platforms directly comparable.
    """

    validation = validate_scale(native_values, contract)

    x, genes = mapped_feature_table(native_values, feature_map)

    if contract.scale in COUNT_SCALES:
        if contract.gene_aggregation != "SUM_BEFORE_TRANSFORM":
            raise ValueError(
                f"{contract.source_expression}: RAW_COUNTS requires "
                "SUM_BEFORE_TRANSFORM"
            )

        gene_counts = aggregate_counts_to_genes(x, genes)
        output = log2_cpm(gene_counts)

        transform_name = "SUM_TO_GENE_THEN_LOG2_CPM_PLUS_1"

    else:
        if contract.gene_aggregation != "MEDIAN_AFTER_TRANSFORM":
            raise ValueError(
                f"{contract.source_expression}: continuous scale requires "
                "MEDIAN_AFTER_TRANSFORM"
            )

        transformed = transform_continuous(x, contract.scale)
        output = aggregate_continuous_to_genes(transformed, genes)

        if contract.scale == ExpressionScale.TPM:
            transform_name = "LOG2_TPM_PLUS_1_THEN_MEDIAN_TO_GENE"

        elif contract.scale == ExpressionScale.FPKM:
            transform_name = "LOG2_FPKM_PLUS_1_THEN_MEDIAN_TO_GENE"

        elif contract.scale in LOG2_OFFSET_SCALES:
            transform_name = "LOG2_X_PLUS_1_THEN_MEDIAN_TO_GENE"

        else:
            transform_name = "IDENTITY_THEN_MEDIAN_TO_GENE"

    report = {
        "source_expression": contract.source_expression,
        "input_scale": contract.scale.value,
        "normalization_provenance": contract.normalization_provenance,
        "native_features": int(native_values.shape[0]),
        "mapped_native_features": int(x.shape[0]),
        "canonical_genes": int(output.shape[0]),
        "samples": int(output.shape[1]),
        "transform": transform_name,
        **validation,
    }

    return output, report
