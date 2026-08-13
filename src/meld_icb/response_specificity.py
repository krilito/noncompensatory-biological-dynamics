"""Response-specific movement estimands."""

from __future__ import annotations

from typing import Any, Mapping, Sequence

from .statistics import auc


def movement_auc(movement_scores: Sequence[float], responder_labels: Sequence[int]) -> float:
    return auc(movement_scores, responder_labels)


def response_specificity_metrics(
    records: Sequence[Mapping[str, Any]],
    movement_rows: Sequence[Mapping[str, Any]],
    *,
    seed: int = 42,
    bootstrap_repeats: int = 2000,
    permutation_repeats: int = 1000,
) -> dict[str, Any]:
    """Run the frozen B3 diagnostic with its small-sample ceiling retained."""
    import numpy as np
    from scipy.stats import mannwhitneyu
    from sklearn.metrics import roc_auc_score

    labels = np.asarray([int(record["y_true"]) for record in records], dtype=int)
    delta_scores = np.asarray([float(row["delta_boundary"]) for row in movement_rows], dtype=float)
    transition_cosines = np.asarray([float(row["cosine_to_reference"]) for row in movement_rows], dtype=float)
    if len(labels) != len(delta_scores):
        raise ValueError("B3 records and movement rows are not aligned")
    if len(set(labels)) < 2:
        return {"analysis_id": "B3_CANONICAL_SPECIFICITY", "canonical_status": "NOT_SUPPORTED", "status_reason": "both response groups are required", "n_records": len(labels)}
    auc_fn = lambda scores, y: float(roc_auc_score(y, scores))
    observed_auc = auc_fn(delta_scores, labels)
    cosine_auc = auc_fn(transition_cosines, labels)
    mw_delta = mannwhitneyu(delta_scores[labels == 1], delta_scores[labels == 0], alternative="two-sided")
    mw_cosine = mannwhitneyu(transition_cosines[labels == 1], transition_cosines[labels == 0], alternative="two-sided")
    rng = np.random.RandomState(seed)
    patient_order: list[Any] = []
    for record in records:
        if record["patient_id"] not in patient_order:
            patient_order.append(record["patient_id"])
    bootstrap: list[float] = []
    for _ in range(bootstrap_repeats):
        drawn = rng.choice(patient_order, size=len(patient_order), replace=True)
        indices = [index for patient in drawn for index, record in enumerate(records) if record["patient_id"] == patient]
        drawn_labels = labels[indices]
        if len(set(drawn_labels.tolist())) == 2:
            bootstrap.append(auc_fn(delta_scores[indices], drawn_labels))
    def permutation_p(scores: np.ndarray) -> float:
        permutation_rng = np.random.RandomState(seed)
        observed = auc_fn(scores, labels)
        count = 0
        for _ in range(permutation_repeats):
            count += int(auc_fn(scores, permutation_rng.permutation(labels)) >= observed)
        return (count + 1) / (permutation_repeats + 1)
    status = "SUPPORTED" if observed_auc > 0.5 and float(mw_delta.pvalue) <= 0.05 and int(np.sum(labels == 0)) >= 2 and int(np.sum(labels == 1)) >= 2 else "NOT_SUPPORTED"
    return {
        "analysis_id": "B3_CANONICAL_SPECIFICITY",
        "n_records": len(labels),
        "n_R": int(np.sum(labels == 1)),
        "n_NR": int(np.sum(labels == 0)),
        "delta_score_auc": float(observed_auc),
        "delta_score_auc_patient_bootstrap_ci": [float(np.percentile(bootstrap, 2.5)), float(np.percentile(bootstrap, 97.5))] if bootstrap else [None, None],
        "delta_score_mannwhitney_p": float(mw_delta.pvalue),
        "delta_score_perm_auc_p": permutation_p(delta_scores),
        "transition_cosine_auc": float(cosine_auc),
        "transition_cosine_mannwhitney_p": float(mw_cosine.pvalue),
        "transition_cosine_perm_auc_p": permutation_p(transition_cosines),
        "canonical_status": status,
        "interpretation": "The small-sample diagnostic did not support response specificity; this is not evidence of biological absence.",
    }
