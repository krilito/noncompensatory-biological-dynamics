"""Paired movement primitives in a fixed representation."""

from __future__ import annotations

from math import sqrt
from typing import Any, Mapping, Sequence


def delta(pre: Mapping[str, float], on_treatment: Mapping[str, float]) -> dict[str, float]:
    if set(pre) != set(on_treatment):
        raise ValueError("paired vectors must have identical axes")
    return {axis: float(on_treatment[axis]) - float(pre[axis]) for axis in pre}


def cosine(left: Mapping[str, float], right: Mapping[str, float]) -> float:
    keys = set(left) & set(right)
    if not keys:
        raise ValueError("vectors have no shared axes")
    dot = sum(float(left[k]) * float(right[k]) for k in keys)
    norm_l = sqrt(sum(float(left[k]) ** 2 for k in keys))
    norm_r = sqrt(sum(float(right[k]) ** 2 for k in keys))
    if norm_l == 0 or norm_r == 0:
        raise ValueError("cosine is undefined for zero vectors")
    return dot / (norm_l * norm_r)


def paired_movement_metrics(
    records: Sequence[Mapping[str, Any]],
    z_axes: Mapping[str, Mapping[str, float]],
    reference_direction: Mapping[str, float],
    freeze: Mapping[str, Any],
) -> dict[str, Any]:
    """Calculate the frozen B2 movement estimand from raw-derived axes."""
    import numpy as np
    from scipy.stats import wilcoxon
    from .frozen_state_transfer import frozen_boundary_score

    rows: list[dict[str, Any]] = []
    for record in records:
        pre_id = str(record["pre_sample_id"])
        edt_id = str(record["edt_sample_id"])
        if pre_id not in z_axes or edt_id not in z_axes:
            raise ValueError(f"paired score mapping is missing {pre_id}/{edt_id}")
        movement = delta(z_axes[pre_id], z_axes[edt_id])
        row = {
            "delta_axes": movement,
            "delta_boundary": float(frozen_boundary_score(z_axes[edt_id], freeze) - frozen_boundary_score(z_axes[pre_id], freeze)),
            "cosine_to_reference": cosine(movement, reference_direction),
        }
        rows.append(row)
    deltas = np.asarray([row["delta_boundary"] for row in rows], dtype=float)
    statistic, p_value = wilcoxon(deltas, alternative="two-sided", zero_method="wilcox")
    cohort_delta = {axis: float(np.mean([row["delta_axes"][axis] for row in rows])) for axis in ("T", "E", "X", "C")}
    fraction = float(np.mean(deltas > 0))
    return {
        "analysis_id": "B2_CANONICAL_MOVEMENT",
        "analysis_unit": "therapy_record",
        "n_records": len(rows),
        "n_patients": len({int(record["patient_id"]) for record in records}),
        "n_delta_gt0": int(np.sum(deltas > 0)),
        "fraction_delta_gt0": fraction,
        "wilcoxon_statistic": float(statistic),
        "wilcoxon_pvalue": float(p_value),
        "cohort_transition_cosine_to_GSE91061": cosine(cohort_delta, reference_direction),
        "median_record_level_cosine": float(np.median([row["cosine_to_reference"] for row in rows])),
        "canonical_status": "SUPPORTED" if len(rows) >= 6 and p_value <= 0.05 and fraction > 0.5 else "NOT_SUPPORTED",
        "interpretation": "movement support does not establish ICB-specific causality",
        "_rows": rows,
    }
