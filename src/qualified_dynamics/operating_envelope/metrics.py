"""Patient aggregation and registered inference for D cohorts."""

from __future__ import annotations

from collections import defaultdict
from typing import Any, Mapping, Sequence

from meld_icb.statistics import registered_auc_inference


class MetricHold(ValueError):
    """Raised when metrics cannot be computed without changing the contract."""


def aggregate_patient_scores(records: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[str, list[Mapping[str, Any]]] = defaultdict(list)
    for record in records:
        patient = str(record.get("patient_id", "")).strip()
        if not patient:
            raise MetricHold("patient aggregation requires patient_id")
        grouped[patient].append(record)
    aggregated: list[dict[str, Any]] = []
    for patient_id, items in sorted(grouped.items()):
        labels = {int(item["response_binary"]) for item in items}
        if len(labels) != 1:
            raise MetricHold(f"patient {patient_id} has conflicting response labels")
        if not items:
            raise MetricHold(f"patient {patient_id} has no eligible biopsies")
        aggregated.append(
            {
                "sample_id": f"patient:{patient_id}",
                "patient_id": patient_id,
                "boundary_score": sum(float(item["boundary_score"]) for item in items) / len(items),
                "response_binary": labels.pop(),
                "biopsy_count": len(items),
            }
        )
    return aggregated


def select_earliest_biopsy(records: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    """Select one biopsy by declared deterministic order only."""

    grouped: dict[str, list[Mapping[str, Any]]] = defaultdict(list)
    for record in records:
        patient = str(record.get("patient_id", "")).strip()
        order = str(record.get("biopsy_order", "")).strip()
        if not patient or not order:
            raise MetricHold("earliest sensitivity requires patient_id and biopsy_order")
        grouped[patient].append(record)
    output: list[dict[str, Any]] = []
    for patient_id, items in sorted(grouped.items()):
        ordered = sorted(items, key=lambda item: (str(item["biopsy_order"]), str(item["sample_id"])))
        first_key = str(ordered[0]["biopsy_order"])
        ties = [item for item in ordered if str(item["biopsy_order"]) == first_key]
        if len(ties) > 1:
            raise MetricHold(f"ambiguous earliest biopsy order for patient {patient_id}")
        output.append(dict(ordered[0]))
    return output


def compute_metrics(records: Sequence[Mapping[str, Any]], *, analysis_unit: str, statistics: Mapping[str, Any]) -> dict[str, Any]:
    if analysis_unit == "patient":
        analysis_records = aggregate_patient_scores(records)
    elif analysis_unit == "sample":
        analysis_records = [dict(record) for record in records]
    else:
        raise MetricHold(f"unsupported analysis unit: {analysis_unit}")
    labels = [int(item["response_binary"]) for item in analysis_records]
    scores = [float(item["boundary_score"]) for item in analysis_records]
    result: dict[str, Any] = {
        "analysis_unit": analysis_unit,
        "n_records": len(analysis_records),
        "n_patients": len({str(item["patient_id"]) for item in records}),
        "labels": {"responder": sum(labels), "nonresponder": len(labels) - sum(labels)},
        "auc": None,
        "ci_95": None,
        "ci_95_low": None,
        "ci_95_high": None,
        "p_two_sided": None,
        "inference_status": "HOLD",
    }
    minimum_class_count = int(statistics.get("minimum_class_count", 1))
    if result["labels"]["responder"] < minimum_class_count or result["labels"]["nonresponder"] < minimum_class_count:
        result["hold_reason"] = "both response classes require the declared minimum count"
        return result
    required_statistics = {
        "inference_contract": "D_OPERATING_ENVELOPE_REGISTERED_V1",
        "bootstrap_resamples": 10000,
        "permutation_resamples": 10000,
        "seed": 20260731,
        "rng_algorithm": "python_random_shared_v1",
        "stream_policy": "shared_bootstrap_then_permutation_persistent_buffer",
        "permutation_comparator": "inclusive_absolute_deviation_ge",
        "correction": "plus_one",
    }
    for key, expected in required_statistics.items():
        if key not in statistics:
            raise MetricHold(f"registered D statistics missing {key}")
        observed = int(statistics[key]) if key in {"bootstrap_resamples", "permutation_resamples", "seed"} else str(statistics[key])
        if observed != expected:
            raise MetricHold(f"registered D statistics require {key}={expected!r}, observed {observed!r}")
    rng_algorithm = str(statistics["rng_algorithm"])
    stream_policy = str(statistics["stream_policy"])
    inference = registered_auc_inference(
        scores,
        labels,
        identifiers=[str(item["sample_id"] if analysis_unit == "sample" else item["patient_id"]) for item in analysis_records],
        identifier_order="sample_id" if analysis_unit == "sample" else "patient_id",
        seed=int(statistics["seed"]),
        bootstrap_resamples=int(statistics["bootstrap_resamples"]),
        permutation_resamples=int(statistics["permutation_resamples"]),
        inference_contract=str(statistics["inference_contract"]),
        rng_algorithm=rng_algorithm,
        stream_policy=stream_policy,
        permutation_comparator=str(statistics["permutation_comparator"]),
        correction=str(statistics["correction"]),
    )
    result.update(
        {
            "auc": float(inference["auc"]),
            "ci_95": tuple(float(value) for value in inference["ci_95"]),
            "ci_95_low": float(inference["ci_95"][0]),
            "ci_95_high": float(inference["ci_95"][1]),
            "p_two_sided": float(inference["p_two_sided"]),
            "bootstrap_resamples": int(inference["bootstrap_resamples"]),
            "permutation_resamples": int(inference["permutation_resamples"]),
            "seed": int(inference["seed"]),
            "rng_algorithm": str(inference["rng_algorithm"]),
            "stream_policy": str(inference["stream_policy"]),
            "inference_contract": str(inference["inference_contract"]),
            "identifier_order": str(inference["identifier_order"]),
            "permutation_comparator": str(inference["permutation_comparator"]),
            "correction": str(inference["correction"]),
            "inference_status": "COMPUTED_REGISTERED",
        }
    )
    return result
