"""Small-sample statistics for the falsification-first paired audit."""

from __future__ import annotations

from itertools import combinations
from typing import Any, Mapping, Sequence

import numpy as np
from scipy import stats

from .statistics import auc


def _as_float(values: Sequence[float]) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    if array.ndim != 1 or not np.isfinite(array).all():
        raise ValueError("scores must be a finite one-dimensional vector")
    return array


def exact_auc_permutation_p(scores: Sequence[float], labels: Sequence[int]) -> float:
    """Enumerate fixed-class label assignments using a two-sided AUC tail."""
    values = _as_float(scores)
    y = np.asarray(labels, dtype=int)
    if values.size != y.size or set(y.tolist()) != {0, 1}:
        raise ValueError("permutation inference requires aligned binary labels")
    observed = auc(values.tolist(), y.tolist())
    n_positive = int(y.sum())
    exceed = 0
    total = 0
    for positive_indices in combinations(range(len(y)), n_positive):
        permuted = np.zeros(len(y), dtype=int)
        permuted[list(positive_indices)] = 1
        candidate = auc(values.tolist(), permuted.tolist())
        exceed += int(abs(candidate - 0.5) >= abs(observed - 0.5) - 1e-15)
        total += 1
    return exceed / total


def _cluster_bootstrap_indices(
    patient_ids: Sequence[str], repeats: int, seed: int
) -> list[list[int]]:
    ordered_patients = list(dict.fromkeys(str(value) for value in patient_ids))
    if len(ordered_patients) < 2:
        raise ValueError("cluster bootstrap requires at least two patients")
    by_patient = {
        patient: [index for index, value in enumerate(patient_ids) if str(value) == patient]
        for patient in ordered_patients
    }
    rng = np.random.default_rng(seed)
    draws: list[list[int]] = []
    for _ in range(repeats):
        selected = rng.choice(ordered_patients, size=len(ordered_patients), replace=True)
        draws.append([index for patient in selected for index in by_patient[str(patient)]])
    return draws


def patient_cluster_bootstrap_indices(
    patient_ids: Sequence[str], repeats: int, seed: int
) -> list[list[int]]:
    """Return deterministic record indices from patient-level resampling."""
    return _cluster_bootstrap_indices(patient_ids, repeats, seed)


def _percentile(values: Sequence[float]) -> list[float | None]:
    finite = np.asarray([value for value in values if np.isfinite(value)], dtype=float)
    if finite.size == 0:
        return [None, None]
    return [float(np.percentile(finite, 2.5)), float(np.percentile(finite, 97.5))]


def _distribution(values: np.ndarray, labels: np.ndarray, label: int) -> dict[str, float | int]:
    group = values[labels == label]
    return {
        "n": int(group.size),
        "median": float(np.median(group)),
        "q1": float(np.percentile(group, 25)),
        "q3": float(np.percentile(group, 75)),
        "mean": float(np.mean(group)),
        "sd": float(np.std(group, ddof=1)) if group.size > 1 else 0.0,
    }


def clustered_auc_analysis(
    channels: Mapping[str, Sequence[float]],
    labels: Sequence[int],
    patient_ids: Sequence[str],
    comparisons: Sequence[tuple[str, str]],
    *,
    bootstrap_repeats: int,
    seed: int,
    compute_record_label_permutation: bool = True,
) -> tuple[dict[str, Any], list[dict[str, float | int]]]:
    """Analyze aligned score channels and their paired AUC differences."""
    arrays = {name: _as_float(values) for name, values in channels.items()}
    y = np.asarray(labels, dtype=int)
    if not arrays or set(y.tolist()) != {0, 1}:
        raise ValueError("clustered AUC analysis requires score channels and both classes")
    if any(len(values) != len(y) for values in arrays.values()) or len(y) != len(patient_ids):
        raise ValueError("clustered AUC inputs must be record-aligned")
    for left, right in comparisons:
        if left not in arrays or right not in arrays:
            raise ValueError(f"unknown AUC comparison: {left} minus {right}")

    observed = {name: float(auc(values.tolist(), y.tolist())) for name, values in arrays.items()}
    channel_draws = {name: [] for name in arrays}
    comparison_draws = {f"{left}_MINUS_{right}": [] for left, right in comparisons}
    draw_rows: list[dict[str, float | int]] = []
    for replicate, indices in enumerate(
        patient_cluster_bootstrap_indices(patient_ids, bootstrap_repeats, seed), start=1
    ):
        drawn_y = y[indices]
        if set(drawn_y.tolist()) != {0, 1}:
            continue
        drawn_auc = {
            name: float(auc(values[indices].tolist(), drawn_y.tolist()))
            for name, values in arrays.items()
        }
        row: dict[str, float | int] = {"replicate": replicate}
        for name, value in drawn_auc.items():
            channel_draws[name].append(value)
            row[f"AUC_{name}"] = value
        for left, right in comparisons:
            key = f"{left}_MINUS_{right}"
            difference = drawn_auc[left] - drawn_auc[right]
            comparison_draws[key].append(difference)
            row[key] = difference
        draw_rows.append(row)

    channel_summary = {}
    for name, values in arrays.items():
        test = stats.mannwhitneyu(values[y == 1], values[y == 0], alternative="two-sided", method="auto")
        channel_summary[name] = {
            "auc": observed[name],
            "auc_ci_95": _percentile(channel_draws[name]),
            "auc_permutation_p_two_sided": (
                exact_auc_permutation_p(values, y) if compute_record_label_permutation else None
            ),
            "mannwhitney_u": float(test.statistic),
            "mannwhitney_p_two_sided": float(test.pvalue),
            "direction_free_auc": max(observed[name], 1.0 - observed[name]),
            "responder": _distribution(values, y, 1),
            "nonresponder": _distribution(values, y, 0),
        }
    comparison_summary = {
        key: {
            "estimate": observed[left] - observed[right],
            "ci_95": _percentile(comparison_draws[key]),
        }
        for (left, right), key in zip(comparisons, comparison_draws)
    }
    return {
        "n_records": len(y),
        "n_patients": len(set(str(value) for value in patient_ids)),
        "n_R": int(y.sum()),
        "n_NR": int(len(y) - y.sum()),
        "bootstrap_repeats_requested": bootstrap_repeats,
        "bootstrap_repeats_valid": len(draw_rows),
        "channels": channel_summary,
        "comparisons": comparison_summary,
        "permutation_unit_note": (
            "exact fixed-class therapy-record label permutation"
            if compute_record_label_permutation
            else "not run: patient-cluster label permutation is undefined because patient 13 has regimen-specific discordant labels"
        ),
    }, draw_rows


def clustered_correlation_analysis(
    x: Sequence[float],
    y: Sequence[float],
    patient_ids: Sequence[str],
    *,
    bootstrap_repeats: int,
    seed: int,
) -> dict[str, Any]:
    """Report Pearson and Spearman correlations with patient-cluster intervals."""
    left = _as_float(x)
    right = _as_float(y)
    if not (len(left) == len(right) == len(patient_ids)):
        raise ValueError("correlation inputs must be record-aligned")
    pearson = stats.pearsonr(left, right)
    spearman = stats.spearmanr(left, right)
    pearson_draws: list[float] = []
    spearman_draws: list[float] = []
    for indices in patient_cluster_bootstrap_indices(patient_ids, bootstrap_repeats, seed):
        if len(np.unique(left[indices])) < 2 or len(np.unique(right[indices])) < 2:
            continue
        pearson_draws.append(float(stats.pearsonr(left[indices], right[indices]).statistic))
        spearman_draws.append(float(stats.spearmanr(left[indices], right[indices]).statistic))
    return {
        "pearson_r": float(pearson.statistic),
        "pearson_p_two_sided": float(pearson.pvalue),
        "pearson_ci_95": _percentile(pearson_draws),
        "spearman_rho": float(spearman.statistic),
        "spearman_p_two_sided": float(spearman.pvalue),
        "spearman_ci_95": _percentile(spearman_draws),
        "bootstrap_repeats_requested": bootstrap_repeats,
        "bootstrap_repeats_valid": min(len(pearson_draws), len(spearman_draws)),
    }


def state_movement_metrics(
    pre_scores: Sequence[float],
    on_scores: Sequence[float],
    labels: Sequence[int],
    patient_ids: Sequence[str],
    *,
    bootstrap_repeats: int,
    seed: int,
) -> dict[str, Any]:
    """Compare PRE, ON, and ON-minus-PRE channels on aligned therapy records."""
    pre = _as_float(pre_scores)
    on = _as_float(on_scores)
    y = np.asarray(labels, dtype=int)
    if not (len(pre) == len(on) == len(y) == len(patient_ids)) or set(y.tolist()) != {0, 1}:
        raise ValueError("state/movement inputs must be aligned and contain both classes")
    movement = on - pre
    channels = {"pre": pre, "on": on, "movement": movement}
    observed_auc = {name: auc(values.tolist(), y.tolist()) for name, values in channels.items()}
    bootstrap = {name: [] for name in channels}
    auc_differences: list[float] = []
    valid = 0
    for indices in _cluster_bootstrap_indices(patient_ids, bootstrap_repeats, seed):
        drawn_y = y[indices]
        if set(drawn_y.tolist()) != {0, 1}:
            continue
        drawn_auc = {
            name: auc(values[indices].tolist(), drawn_y.tolist())
            for name, values in channels.items()
        }
        for name, value in drawn_auc.items():
            bootstrap[name].append(value)
        auc_differences.append(drawn_auc["on"] - drawn_auc["movement"])
        valid += 1
    result: dict[str, Any] = {
        "n_records": len(y),
        "n_patients": len(set(str(value) for value in patient_ids)),
        "n_R": int(y.sum()),
        "n_NR": int(len(y) - y.sum()),
        "bootstrap_repeats_requested": bootstrap_repeats,
        "bootstrap_repeats_valid": valid,
        "on_minus_movement_auc": float(observed_auc["on"] - observed_auc["movement"]),
        "on_minus_movement_auc_ci_95": _percentile(auc_differences),
    }
    for name, values in channels.items():
        responder = values[y == 1]
        nonresponder = values[y == 0]
        test = stats.mannwhitneyu(responder, nonresponder, alternative="two-sided", method="auto")
        result[name] = {
            "auc": float(observed_auc[name]),
            "auc_ci_95": _percentile(bootstrap[name]),
            "auc_permutation_p_two_sided": exact_auc_permutation_p(values, y),
            "mannwhitney_u": float(test.statistic),
            "mannwhitney_p_two_sided": float(test.pvalue),
            "responder": _distribution(values, y, 1),
            "nonresponder": _distribution(values, y, 0),
        }
    return result


def benjamini_hochberg(p_values: Sequence[float]) -> list[float]:
    values = np.asarray(p_values, dtype=float)
    order = np.argsort(values)
    adjusted = np.empty(len(values), dtype=float)
    running = 1.0
    for reverse_rank, index in enumerate(order[::-1], start=1):
        rank = len(values) - reverse_rank + 1
        running = min(running, float(values[index]) * len(values) / rank)
        adjusted[index] = running
    return adjusted.tolist()


def paired_change_control(
    pre_scores: Sequence[float],
    on_scores: Sequence[float],
    patient_ids: Sequence[str],
    *,
    bootstrap_repeats: int,
    seed: int,
) -> dict[str, Any]:
    pre = _as_float(pre_scores)
    on = _as_float(on_scores)
    if not (len(pre) == len(on) == len(patient_ids)):
        raise ValueError("paired control inputs must be aligned")
    change = on - pre
    test = stats.wilcoxon(change, alternative="two-sided", zero_method="wilcox")
    ranks = stats.rankdata(np.abs(change[change != 0]))
    nonzero = change[change != 0]
    positive_rank = float(ranks[nonzero > 0].sum())
    negative_rank = float(ranks[nonzero < 0].sum())
    denominator = positive_rank + negative_rank
    boot_medians = [
        float(np.median(change[indices]))
        for indices in _cluster_bootstrap_indices(patient_ids, bootstrap_repeats, seed)
    ]
    return {
        "n_records": len(change),
        "n_patients": len(set(str(value) for value in patient_ids)),
        "median_change": float(np.median(change)),
        "mean_change": float(np.mean(change)),
        "median_change_ci_95": _percentile(boot_medians),
        "wilcoxon_statistic": float(test.statistic),
        "wilcoxon_p_two_sided": float(test.pvalue),
        "paired_rank_biserial": (positive_rank - negative_rank) / denominator if denominator else 0.0,
        "n_positive": int(np.sum(change > 0)),
        "n_negative": int(np.sum(change < 0)),
        "n_zero": int(np.sum(change == 0)),
        "fraction_positive": float(np.mean(change > 0)),
    }


def change_coupling(pre_scores: Sequence[float], on_scores: Sequence[float]) -> dict[str, float]:
    pre = _as_float(pre_scores)
    on = _as_float(on_scores)
    movement = on - pre
    average = (on + pre) / 2.0
    pearson_pre = stats.pearsonr(pre, movement)
    spearman_pre = stats.spearmanr(pre, movement)
    spearman_average = stats.spearmanr(average, movement)
    return {
        "pearson_pre_vs_delta_r": float(pearson_pre.statistic),
        "pearson_pre_vs_delta_p": float(pearson_pre.pvalue),
        "spearman_pre_vs_delta_rho": float(spearman_pre.statistic),
        "spearman_pre_vs_delta_p": float(spearman_pre.pvalue),
        "oldham_average_vs_delta_rho": float(spearman_average.statistic),
        "oldham_average_vs_delta_p": float(spearman_average.pvalue),
    }
