"""Machine-readable contracts for the D operating-envelope lane."""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path
from typing import Any, Mapping


class ContractError(ValueError):
    """Raised when a cohort contract cannot be safely executed."""


VALID_ANALYSIS_UNITS = {"sample", "patient"}
VALID_REPRESENTATIONS = {"gene_expression", "canonical_axis_z"}
VALID_STATUSES = {"AUTO", "SUPPORTED", "NEAR_NULL", "INVERTED", "UNRESOLVED", "HOLD"}
INFERENCE_CONTRACTS = {
    "D_OPERATING_ENVELOPE_REGISTERED_V1": {
        "bootstrap_resamples": 10000,
        "permutation_resamples": 10000,
        "seed": 20260731,
        "rng_algorithm": "python_random_shared_v1",
        "stream_policy": "shared_bootstrap_then_permutation_persistent_buffer",
    }
}
REFERENCE_ONLY_MARKERS = (
    "final_numeric_registry",
    "canonical_result_registry",
    "panel_source_data",
    "producer_outputs",
    "source_data",
    "figure_source",
)

REQUIRED_KEYS = {
    "dataset_id",
    "dataset_contract",
    "local_config_key",
    "input_files",
    "input_adapter",
    "inclusion_rules",
    "exclusion_rules",
    "analysis_unit",
    "response",
    "scorer",
    "statistics",
    "carrier_qualification",
    "verdict",
    "provenance",
}


def _load_json_yaml(path: Path) -> dict[str, Any]:
    """Load JSON-compatible YAML without adding a runtime dependency."""
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ContractError(f"unreadable contract: {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise ContractError(f"contract must be an object: {path}")
    return value


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _contains_reference_marker(value: Any) -> bool:
    if isinstance(value, Mapping):
        return any(_contains_reference_marker(k) or _contains_reference_marker(v) for k, v in value.items())
    if isinstance(value, (list, tuple, set)):
        return any(_contains_reference_marker(item) for item in value)
    text = str(value).lower().replace(chr(92), "/")
    return any(marker in text for marker in REFERENCE_ONLY_MARKERS)


def validate_contract(contract: Mapping[str, Any]) -> dict[str, Any]:
    """Validate and normalize one D contract.

    Validation is deliberately conservative: every executable contract must
    carry input/provenance declarations and a declared verdict.  A contract
    may still be ``HOLD``; that is an explicit scientific state, not an error.
    """
    missing = sorted(REQUIRED_KEYS.difference(contract))
    if missing:
        raise ContractError(f"contract missing keys: {missing}")
    dataset_id = str(contract["dataset_id"]).strip()
    if not dataset_id:
        raise ContractError("dataset_id must be non-empty")
    unit = str(contract["analysis_unit"]).strip().lower()
    if unit not in VALID_ANALYSIS_UNITS:
        raise ContractError(f"unsupported analysis unit: {unit}")
    files = contract["input_files"]
    if not isinstance(files, list) or not files or any(not str(item).strip() for item in files):
        raise ContractError("input_files must be a non-empty list")
    for item in files:
        if Path(str(item)).is_absolute():
            raise ContractError("input_files must be relative/local-config paths")
    adapter = str(contract["input_adapter"]).strip()
    if not adapter:
        raise ContractError("input_adapter must be declared")
    response = contract["response"]
    if not isinstance(response, Mapping):
        raise ContractError("response must be an object")
    labels = response.get("positive_labels"), response.get("negative_labels")
    if not all(isinstance(item, list) and item for item in labels):
        raise ContractError("response positive_labels and negative_labels are required")
    if set(map(str.upper, labels[0])) & set(map(str.upper, labels[1])):
        raise ContractError("response classes overlap")
    scorer = contract["scorer"]
    if not isinstance(scorer, Mapping) or scorer.get("authority") != "meld_icb.frozen_state_transfer":
        raise ContractError("D must reuse the single A scorer authority")
    if scorer.get("preprocessing") != "PRE_MEAN_THEN_AXIS_Z" or int(scorer.get("ddof", 0)) != 1:
        raise ContractError("D scorer preprocessing/ddof diverges from frozen A")
    stats = contract["statistics"]
    if not isinstance(stats, Mapping):
        raise ContractError("statistics must be an object")
    for key in ("seed", "bootstrap_resamples", "permutation_resamples", "ci", "auc", "rng_algorithm", "stream_policy", "inference_contract", "permutation_comparator", "correction"):
        if key not in stats:
            raise ContractError(f"statistics missing {key}")
    inference_contract = str(stats["inference_contract"]).strip()
    expected_stats = INFERENCE_CONTRACTS.get(inference_contract)
    if expected_stats is None:
        raise ContractError(f"unsupported inference contract: {inference_contract}")
    for key, expected in expected_stats.items():
        observed = int(stats[key]) if key in {"seed", "bootstrap_resamples", "permutation_resamples"} else str(stats[key])
        if observed != expected:
            raise ContractError(f"D inference contract {inference_contract} requires {key}={expected!r}")
    if str(stats["permutation_comparator"]) != "inclusive_absolute_deviation_ge":
        raise ContractError("D requires inclusive_absolute_deviation_ge comparator")
    if str(stats["correction"]) != "plus_one":
        raise ContractError("D requires plus_one correction")
    verdict = contract["verdict"]
    if not isinstance(verdict, Mapping):
        raise ContractError("verdict must be an object")
    declared = str(verdict.get("declared_status", "")).upper()
    if declared not in VALID_STATUSES:
        raise ContractError(f"invalid declared verdict: {declared}")
    if bool(verdict.get("posthoc_inversion", False)):
        raise ContractError("post hoc score inversion is prohibited")
    inclusion = contract["inclusion_rules"]
    exclusion = contract["exclusion_rules"]
    if not isinstance(inclusion, Mapping) or not isinstance(exclusion, Mapping):
        raise ContractError("inclusion_rules and exclusion_rules must be objects")
    selection = str(inclusion.get("biopsy_selection", "all_eligible")).lower()
    if selection not in {"all_eligible", "earliest_sensitivity", "not_applicable"}:
        raise ContractError("biopsies cannot be selected by response, score, threshold, or outcome")
    if any(token in selection for token in ("score", "response", "outcome", "threshold")):
        raise ContractError("outcome/score-based biopsy selection is prohibited")
    provenance = contract["provenance"]
    if not isinstance(provenance, Mapping):
        raise ContractError("provenance must be an object")
    source_path = str(provenance.get("lineage_source_path", "")).strip()
    source_hash = str(provenance.get("lineage_source_sha256", "")).strip().lower()
    if not source_path or not re.fullmatch(r"[0-9a-f]{64}", source_hash):
        raise ContractError("complete lineage source path and SHA256 are required")
    source_path_obj = Path(source_path)
    if source_path_obj.is_absolute() or ".." in source_path_obj.parts:
        raise ContractError("lineage_source_path must be a relative canonical authority path")
    if _contains_reference_marker(source_path):
        raise ContractError("reference-only objects cannot be producer authority")
    normalized = dict(contract)
    normalized.update({"dataset_id": dataset_id, "analysis_unit": unit})
    normalized["verdict"] = {**dict(verdict), "declared_status": declared, "posthoc_inversion": False}
    return normalized


def load_contract(path: str | Path) -> dict[str, Any]:
    contract = _load_json_yaml(Path(path))
    return validate_contract(contract)


def contract_sha256(path: str | Path) -> str:
    return _sha256(Path(path))
