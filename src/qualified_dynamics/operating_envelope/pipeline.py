"""The sole D operating-envelope execution pipeline."""

from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
from typing import Any, Mapping

from meld_icb.frozen_state_transfer import (
    AXIS_ORDER,
    frozen_boundary_score,
    pre_mean_then_axis_z,
)

from .contracts import (
    REFERENCE_ONLY_MARKERS,
    ContractError,
    contract_sha256,
    load_contract,
)
from .inputs import InputHold, encode_response, load_rows, validate_rows
from .metrics import MetricHold, compute_metrics, select_earliest_biopsy
from .verdicts import classify_verdict


def _logical_path(path: Path, root: Path) -> str:
    """Render a repository-relative/logical path without local absolute roots."""
    try:
        return path.resolve().relative_to(root.resolve()).as_posix()
    except ValueError:
        return path.name


def _sha256(path: Path) -> str:
    return contract_sha256(path)


def _write(payload: Mapping[str, Any], destination: Path) -> dict[str, Any]:
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n", encoding="utf-8")
    return dict(payload)


def _source_binding(contract: Mapping[str, Any], source_root: Path | None) -> dict[str, Any]:
    provenance = contract["provenance"]
    declared = str(provenance["lineage_source_path"])
    expected = str(provenance["lineage_source_sha256"]).lower()
    if source_root is None:
        return {
            "status": "UNVERIFIED_LINEAGE",
            "reason": "lineage source was not locally cross-checked; public execution remains valid",
            "declared_path": declared,
            "expected_sha256": expected,
        }
    source = source_root / declared
    if not source.exists() or not source.is_file():
        return {"status": "HOLD", "reason": "lineage source is missing", "declared_path": declared, "expected_sha256": expected}
    observed = _sha256(source)
    if observed != expected:
        return {"status": "HOLD", "reason": "lineage source SHA256 mismatch", "declared_path": declared, "expected_sha256": expected, "observed_sha256": observed}
    return {"status": "PASS", "declared_path": declared, "expected_sha256": expected, "observed_sha256": observed}


def _resolve_input(input_path: str | Path | None, contract: Mapping[str, Any], data_root: Path | None) -> Path:
    declared = Path(str(contract["input_files"][0]))
    if declared.is_absolute():
        raise InputHold("contract input path must remain relative")
    if input_path is not None:
        candidate = Path(input_path)
        if candidate.is_dir():
            return candidate / declared
        return candidate
    root = data_root or Path(os.environ.get("MELD_DATA_ROOT", "data/external"))
    return root / declared


def _reference_only(path: Path) -> bool:
    lowered = path.as_posix().lower()
    return any(marker in lowered for marker in REFERENCE_ONLY_MARKERS)


def _transform_expression(row: Mapping[str, Any], contract: Mapping[str, Any]) -> dict[str, float]:
    values = row.get("_expression", {})
    if not isinstance(values, Mapping) or not values:
        raise InputHold(f"sample {row.get('sample_id', '<unknown>')} has no expression")
    transform = str(contract["scorer"].get("expression_transform", "")).lower()
    if transform != "truncate_negative_then_log2p1":
        raise InputHold(f"unsupported declared expression transform: {transform}")
    if row.get("_expression_transformed"):
        return {str(key): float(value) for key, value in values.items()}
    return {str(key): math.log2(max(0.0, float(value)) + 1.0) for key, value in values.items()}


def _score_rows(rows: list[dict[str, Any]], contract: Mapping[str, Any], freeze: Mapping[str, Any], representation: str) -> list[dict[str, Any]]:
    if representation == "canonical_axis_z":
        axes = [{axis: float(row[axis]) for axis in AXIS_ORDER} for row in rows]
    else:
        expression = [_transform_expression(row, contract) for row in rows]
        axes, _, _ = pre_mean_then_axis_z(
            expression,
            freeze,
            [str(row["sample_id"]) for row in rows],
            scope_id=f"D:{contract['dataset_id']}:cohort",
        )
    scored: list[dict[str, Any]] = []
    for row, axis in zip(rows, axes):
        scored.append(
            {
                "sample_id": str(row["sample_id"]),
                "patient_id": str(row["patient_id"]),
                "timepoint": str(row.get("timepoint", "")),
                "biopsy_order": str(row.get("biopsy_order", "")),
                **axis,
                "boundary_score": frozen_boundary_score(axis, freeze),
                "response_binary": encode_response(row[contract["response"]["column"]], contract),
            }
        )
    return scored


def _summary(rows: list[Mapping[str, Any]]) -> dict[str, Any]:
    patients: dict[str, list[str]] = {}
    for row in rows:
        patients.setdefault(str(row["patient_id"]), []).append(str(row["sample_id"]))
    return {
        "n_biopsies": len(rows),
        "n_patients": len(patients),
        "n_duplicate_patients": sum(len(samples) > 1 for samples in patients.values()),
    }


def _compare_counts(summary: Mapping[str, Any], expected: Mapping[str, Any]) -> dict[str, Any]:
    checks: dict[str, bool] = {}
    if "n_biopsies" in expected:
        checks["n_biopsies"] = int(summary["n_biopsies"]) == int(expected["n_biopsies"])
    if "n_patients" in expected:
        checks["n_patients"] = int(summary["n_patients"]) == int(expected["n_patients"])
    if "n_duplicate_patients" in expected:
        checks["n_duplicate_patients"] = int(summary["n_duplicate_patients"]) == int(expected["n_duplicate_patients"])
    if not checks:
        return {"status": "NOT_DECLARED", "checks": checks}
    return {"status": "PASS" if all(checks.values()) else "MISMATCH", "checks": checks}


def run_operating_envelope(
    *,
    config_path: str | Path,
    input_path: str | Path | None = None,
    output_path: str | Path,
    representation: str | None = None,
    analysis_unit: str | None = None,
    source_root: str | Path | None = None,
    repo_root: str | Path | None = None,
    data_root: str | Path | None = None,
) -> dict[str, Any]:
    config_file = Path(config_path)
    destination = Path(output_path)
    try:
        contract = load_contract(config_file)
    except (ContractError, OSError) as exc:
        return _write({"status": "HOLD", "verdict": "HOLD", "hold_reason": str(exc)}, destination)
    source_binding = _source_binding(contract, Path(source_root) if source_root is not None else None)
    root = Path(repo_root) if repo_root is not None else config_file.parents[2]
    source = _resolve_input(input_path, contract, Path(data_root) if data_root is not None else None)
    authority_path = root / "src" / "qualified_dynamics" / "operating_envelope" / "pipeline.py"
    metadata: dict[str, Any] = {
        "config_path": _logical_path(config_file, root),
        "config_sha256": _sha256(config_file) if config_file.exists() else None,
        "input_path": _logical_path(source, Path(data_root) if data_root is not None else root),
        "input_representation": representation or contract["scorer"].get("representation", "gene_expression"),
        "analysis_unit": analysis_unit or contract["analysis_unit"],
        "lineage_source": source_binding,
        "public_authority": {
            "module": "qualified_dynamics.operating_envelope.pipeline",
            "path": _logical_path(authority_path, root),
            "sha256": _sha256(authority_path) if authority_path.exists() else None,
            "status": "PASS" if authority_path.exists() else "HOLD",
        },
        "authority": "qualified_dynamics.operating_envelope.pipeline",
    }
    if source_binding.get("status") == "HOLD":
        return _write({"status": "HOLD", "verdict": "HOLD", "hold_reason": source_binding.get("reason", "producer source unverified"), "metadata": metadata}, destination)
    if _reference_only(source):
        return _write({"status": "HOLD", "verdict": "HOLD", "hold_reason": "reference-only object cannot be producer input", "metadata": metadata}, destination)
    if not source.exists():
        return _write({"status": "HOLD", "verdict": "HOLD", "hold_reason": "required cohort input file is missing", "metadata": metadata}, destination)
    if str(contract["verdict"]["declared_status"]).upper() == "HOLD":
        return _write({"status": "HOLD", "verdict": "HOLD", "hold_reason": contract.get("hold_reason", contract["verdict"].get("primary_rule", "cohort contract is explicitly held")), "metadata": metadata}, destination)
    unit = analysis_unit or str(contract["analysis_unit"])
    rep = representation or str(contract["scorer"].get("representation", "gene_expression"))
    try:
        rows = load_rows(input_path=source, contract=contract, data_root=Path(data_root) if data_root is not None else None)
        validate_rows(rows, contract, representation=rep)
        freeze_path = root / "configs" / "frozen_state_transfer.yaml"
        from meld_icb.config import load_freeze
        freeze = load_freeze(root)
        scored = _score_rows(rows, contract, freeze, rep)
        metrics = compute_metrics(scored, analysis_unit=unit, statistics=contract["statistics"])
        populations: dict[str, Any] = {"primary": {"metrics": metrics, "scored_rows": scored}}
        if unit == "patient" and len({str(row["patient_id"]) for row in scored}) < len(scored):
            try:
                earliest = select_earliest_biopsy(scored)
                populations["earliest_sensitivity"] = {"metrics": compute_metrics(earliest, analysis_unit="patient", statistics=contract["statistics"]), "scored_rows": earliest}
            except MetricHold as exc:
                populations["earliest_sensitivity"] = {"status": "HOLD", "hold_reason": str(exc)}
            populations["biopsy_clustered_sensitivity"] = {"status": "HOLD", "hold_reason": "clustered biopsy inference is sensitivity-only and requires a registered cluster estimator"}
        summary = _summary(scored)
        comparison = _compare_counts(summary, contract.get("expected_counts", {}))
        verdict = classify_verdict(metrics, contract["verdict"])
        payload = {
            "status": "COMPUTED" if metrics["inference_status"] == "COMPUTED_REGISTERED" else "HOLD",
            "verdict": verdict if metrics["inference_status"] == "COMPUTED_REGISTERED" else "HOLD",
            "cohort_id": contract["dataset_id"],
            "populations": populations,
            "population_summary": summary,
            "reference_comparison": comparison,
            "orientation": {"declared": contract["scorer"]["orientation"], "posthoc_inversion": False},
            "metadata": {**metadata, "input_sha256": _sha256(source), "preprocessing_id": contract["scorer"]["preprocessing"], "workbook_sheet": contract.get("dataset_contract", {}).get("workbook_sheet")},
        }
        if comparison.get("status") == "MISMATCH":
            payload["status"] = "HOLD"
            payload["verdict"] = "HOLD"
            payload["hold_reason"] = "declared population count contract mismatch"
    except (InputHold, MetricHold, ContractError, ValueError, OSError, KeyError) as exc:
        payload = {"status": "HOLD", "verdict": "HOLD", "hold_reason": str(exc), "metadata": metadata}
    return _write(payload, destination)


def run_declared_cohorts(
    *,
    config_dir: str | Path,
    data_root: str | Path,
    output_dir: str | Path,
    source_root: str | Path | None = None,
    repo_root: str | Path | None = None,
    representation: str | None = None,
) -> dict[str, Any]:
    outputs: dict[str, Any] = {}
    config_paths = sorted(Path(config_dir).glob("*.yaml"))
    for config_path in config_paths:
        contract = load_contract(config_path)
        outputs[contract["dataset_id"]] = run_operating_envelope(
            config_path=config_path,
            input_path=Path(data_root),
            output_path=Path(output_dir) / f"{contract['dataset_id']}.json",
            representation=representation,
            analysis_unit=str(contract["analysis_unit"]),
            source_root=source_root,
            repo_root=repo_root,
            data_root=data_root,
        )
    overall = "COMPUTED" if outputs and all(item.get("status") == "COMPUTED" for item in outputs.values()) else "HOLD"
    return {"producer": "D_OPERATING_ENVELOPE", "status": overall, "authority": "qualified_dynamics.operating_envelope.pipeline", "cohorts": outputs, "data_root": _logical_path(Path(data_root), Path(repo_root) if repo_root is not None else Path(config_dir).parents[2])}
