"""Public figure/table adapters over the four canonical shared chains.

Adapters deliberately perform only projection, ordering, formatting, and
expected-value validation.  Scientific statistics are produced once by the
chain runners; this module never imports a statistics routine or raw cohort
loader.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Mapping

def load_object_map(repo_root: str | Path) -> dict[str, dict[str, Any]]:
    """Load and validate the independent object-to-chain contract."""
    path = Path(repo_root) / "manifests" / "object_adapter_map.json"
    payload = _read_json(path)
    if not payload or payload.get("schema_version") != "OBJECT_ADAPTER_MAP_V1" or not isinstance(payload.get("objects"), Mapping):
        raise ValueError(f"invalid object adapter map: {path}")
    objects: dict[str, dict[str, Any]] = {}
    for object_id, raw in payload["objects"].items():
        if not isinstance(raw, Mapping) or raw.get("chain") not in {"A", "B", "C", "D"}:
            raise ValueError(f"invalid adapter entry: {object_id}")
        disposition = str(raw.get("disposition", ""))
        default_status = str(raw.get("default_status", ""))
        kind = (
            "CONCEPTUAL" if disposition == "CONCEPTUAL_ASSET_CONTRACT"
            else "EXCLUDED" if disposition == "EXCLUDED_BY_POLICY"
            else "RETIRED" if disposition == "RETIRED" or default_status in {"RETIRED", "SUPERSEDED_HOLD"}
            else "MIXED" if default_status == "COMPUTED_PRIMARY_PLUS_COMPARISON_ONLY"
            else "MAPPED" if default_status == "COMPUTED"
            else "HOLD"
        )
        objects[str(object_id)] = {**dict(raw), "kind": kind, "authority": default_status}
    if len(objects) != 25:
        raise ValueError(f"expected 25 object adapters, observed {len(objects)}")
    return objects

CHAIN_FILES = {
    "A": ("A_PUBLIC_NUMERIC_OUTPUT.json", "A_85_NUMERIC_OUTPUT.json"),
    "B": ("B_PUBLIC_NUMERIC_OUTPUT.json",),
    "C": ("C_PUBLIC_NUMERIC_OUTPUT.json",),
    "D": ("D_PUBLIC_NUMERIC_OUTPUT.json", "D_OPERATING_ENVELOPE_OUTPUT.json"),
}


def _read_json(path: Path) -> dict[str, Any] | None:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    return value if isinstance(value, dict) else None


def _chain_payload(shared_dir: Path, chain: str) -> dict[str, Any] | None:
    for name in CHAIN_FILES[chain]:
        payload = _read_json(shared_dir / name)
        if payload is not None:
            return payload
    return None


def _status_from_chain(payload: Mapping[str, Any] | None) -> str:
    if payload is None:
        return "HOLD"
    status = str(payload.get("status", payload.get("verdict", "HOLD"))).upper()
    return "COMPUTED" if status in {"COMPUTED", "GATES_PASSED"} else "HOLD"


def _b_projection(chain: Mapping[str, Any]) -> dict[str, Any]:
    # Explicit allow-list prevents row-level records and source IDs reaching
    # figure/table source data.
    b2 = chain.get("b2_movement", {})
    b3 = chain.get("b3_response_specificity", {})
    b4 = chain.get("b4_incremental_foldsafe", {})
    keys = {
        "b2": ("n_records", "n_patients", "n_delta_gt0", "fraction_delta_gt0", "median_record_level_cosine", "cohort_transition_cosine_to_GSE91061", "wilcoxon_pvalue", "canonical_status"),
        "b3": ("n_records", "n_R", "n_NR", "delta_score_auc", "delta_score_auc_patient_bootstrap_ci", "delta_score_mannwhitney_p", "delta_score_perm_auc_p", "transition_cosine_auc", "transition_cosine_mannwhitney_p", "transition_cosine_perm_auc_p", "canonical_status"),
        "b4": ("n_records", "n_patients", "delta_L", "L_baseline_X0_pre", "L_augmented_X1_pre_plus_delta", "P_improvement", "P_worsening", "P_two_sided", "n_sign_patterns", "oof_auc_X0", "oof_auc_X1", "secondary_delta_auc", "canonical_status"),
    }
    return {name: {key: source[key] for key in fields if key in source} for name, source, fields in (("b2", b2, keys["b2"]), ("b3", b3, keys["b3"]), ("b4", b4, keys["b4"]))}


def _c_projection(chain: Mapping[str, Any]) -> dict[str, Any]:
    gates = chain.get("gates", {})
    gate_values = {
        "n_cells": chain.get("n_cells", gates.get("n_on_treatment_cells", {}).get("observed")),
        "n_signatures": chain.get("n_signatures"),
        "lymphoid_share": gates.get("A_lymphoid_share", {}).get("observed"),
        "malignant_share": gates.get("B_malignant_share", {}).get("observed"),
        "top_two_carriers": gates.get("C_top_two_carriers", {}).get("observed"),
    }
    summary = []
    for row in chain.get("summary", []):
        if not isinstance(row, Mapping):
            continue
        summary.append({key: row[key] for key in ("signature", "intended_referent", "dominant_carrier", "lymphoid_share", "myeloid_share", "immune_share", "malignant_share", "intent_carrier_concordant", "discordant_tumour_intrinsic") if key in row})
    gate_values["signature_summary"] = summary
    return gate_values


def _aggregate_projection(chain: Mapping[str, Any]) -> dict[str, Any]:
    if "cohorts" in chain and isinstance(chain["cohorts"], Mapping):
        rows = []
        for cohort, item in sorted(chain["cohorts"].items()):
            if not isinstance(item, Mapping):
                continue
            row: dict[str, Any] = {"cohort": str(cohort), "status": str(item.get("status", "HOLD"))}
            metrics = item.get("metrics", {})
            if not isinstance(metrics, Mapping):
                metrics = {}
            if not metrics and isinstance(item.get("populations"), Mapping):
                primary = item["populations"].get("primary", {})
                if isinstance(primary, Mapping) and isinstance(primary.get("metrics"), Mapping):
                    metrics = primary["metrics"]
            if isinstance(metrics, Mapping):
                for key in ("auc", "ci_95", "p_two_sided", "n_records", "n_patients", "inference_status"):
                    if key in metrics:
                        row[key] = metrics[key]
            for key in ("verdict", "population_summary", "reference_comparison"):
                if key in item and key == "verdict":
                    row[key] = item[key]
            rows.append(row)
        projection: dict[str, Any] = {"rows": rows}
        pooled = chain.get("pooled")
        if isinstance(pooled, Mapping):
            projection["pooled"] = {key: pooled[key] for key in ("analysis_unit", "n_records", "auc", "ci_95", "p_two_sided", "inference_status") if key in pooled}
        identity = chain.get("identity")
        if isinstance(identity, Mapping):
            projection["identity"] = {key: identity[key] for key in ("n_records", "sha256") if key in identity}
        return projection
    return {key: chain[key] for key in ("status", "verdict", "auc", "ci_95", "p_two_sided", "n_records", "n_patients") if key in chain}


def _expected_validation(object_id: str, aggregate: Mapping[str, Any]) -> dict[str, Any]:
    checks: dict[str, bool] = {}
    details: dict[str, dict[str, Any]] = {}
    def check(name: str, observed: Any, expected: Any, tolerance: float = 1e-12) -> None:
        try:
            error = abs(float(observed) - float(expected))
            passed = math.isfinite(error) and error <= tolerance
        except (TypeError, ValueError):
            error, passed = None, False
        checks[name] = passed
        details[name] = {"observed": observed, "expected": expected, "error": error, "tolerance": tolerance}
    if object_id in {"Figure2", "S05", "S06", "S07", "S08", "S11", "ST02"}:
        rows = {str(row.get("cohort")): row for row in aggregate.get("rows", []) if isinstance(row, Mapping)}
        check("PRJEB23709_auc", rows.get("PRJEB23709", {}).get("auc"), 0.8441558441558441)
        check("MORRISON_auc", rows.get("MORRISON-1-public", {}).get("auc"), 0.825)
        check("MGH_auc", rows.get("MGH", {}).get("auc"), 0.675)
        pooled = aggregate.get("pooled", {})
        identity = aggregate.get("identity", {})
        check("pooled_n_records", pooled.get("n_records"), 85)
        check("pooled_auc", pooled.get("auc"), 0.7705345501955672)
        check("pooled_p_two_sided", pooled.get("p_two_sided"), 0.000999000999000999)
        ci = pooled.get("ci_95") or [None, None]
        check("pooled_ci_low", ci[0], 0.6538461538461539)
        check("pooled_ci_high", ci[1], 0.8800521512385919)
        check("identity_n_records", identity.get("n_records"), 85)
        identity_expected = "bedb09abc0ee1338fc904f834cda957751a79202cff9102f1c5fe986a11ddd11"
        identity_observed = identity.get("sha256")
        identity_passed = identity_observed == identity_expected
        checks["identity_sha256"] = identity_passed
        details["identity_sha256"] = {"observed": identity_observed, "expected": identity_expected, "error": 0 if identity_passed else 1, "comparison_mode": "exact_string", "tolerance": 0}
    elif object_id in {"Figure5", "S09", "S10", "ST06"}:
        rows = {str(row.get("cohort")): row for row in aggregate.get("rows", []) if isinstance(row, Mapping)}
        gse = rows.get("GSE78220", {})
        check("GSE78220_auc", gse.get("auc"), 0.47619047619047616)
        check("GSE78220_ci_low", (gse.get("ci_95") or [None, None])[0], 0.24404761904761904)
        check("GSE78220_ci_high", (gse.get("ci_95") or [None, None])[1], 0.7083333333333334)
        check("GSE78220_p", gse.get("p_two_sided"), 0.8573142685731426)
        check("GSE78220_n_patients", gse.get("n_patients"), 26)
    elif object_id in {"Figure3", "S03", "ST03"}:
        b2, b3, b4 = aggregate.get("b2", {}), aggregate.get("b3", {}), aggregate.get("b4", {})
        check("b2_n_records", b2.get("n_records"), 16)
        check("b2_n_patients", b2.get("n_patients"), 15)
        check("b2_median_cosine", b2.get("median_record_level_cosine"), 0.6993631612147504)
        check("b3_delta_auc", b3.get("delta_score_auc"), 0.3272727272727272)
        check("b4_two_sided_p", b4.get("P_two_sided"), 0.03094482421875)
    elif object_id in {"Figure4", "S12", "ST04"}:
        check("n_cells", aggregate.get("n_cells"), 10363)
        check("lymphoid_share", aggregate.get("lymphoid_share"), 0.8629283905029297)
        check("malignant_share", aggregate.get("malignant_share"), 0.0007919137715362012)
    if details and not any(detail.get("observed") is not None for detail in details.values()):
        return {"status": "NOT_APPLICABLE", "checks": {}, "details": {}, "reason": "required aggregate fields are absent"}
    return {"status": "PASS" if checks and all(checks.values()) else "MISMATCH", "checks": checks, "details": details, "tolerance": 1e-12}


def build_adapter_payload(object_id: str, *, shared_dir: str | Path, repo_root: str | Path) -> dict[str, Any]:
    object_map = load_object_map(repo_root)
    if object_id not in object_map:
        raise KeyError(f"unknown object_id: {object_id}")
    spec = object_map[object_id]
    shared = Path(shared_dir)
    root = Path(repo_root)
    base: dict[str, Any] = {"schema_version": "OBJECT_ADAPTER_V1", "object_id": object_id, "shared_chain": spec["chain"], "authority_status": spec["authority"], "adapter_only": True}
    if spec["kind"] == "CONCEPTUAL":
        base.update({"status": "CONCEPTUAL_ASSET_CONTRACT", "conceptual_asset": {"type": "F1_PANEL_CONTRACT", "source_rows": "aggregate panel rows only", "cohort_roles": ["derivation", "external_transfer", "carrier_context"], "questions": ["frozen coefficients", "state transfer", "paired movement", "response specificity", "incremental value"]}})
        return base
    if spec["kind"] == "EXCLUDED":
        base.update({"status": "EXCLUDED_BY_POLICY", "exclusion_reason": "identifier/license mapping is not closed", "schema_only": {"omitted_patient_details_included": False, "internal_ids_included": False, "hmac_mapping": "not available"}})
        return base
    if spec["kind"] == "RETIRED":
        base.update({"status": "RETIRED", "retirement_reason": spec.get("retirement_reason", "superseded object is excluded from active computation"), "active_projection": False})
        return base
    chain = _chain_payload(shared, spec["chain"])
    chain_status = _status_from_chain(chain)
    base["chain_status"] = chain_status
    if chain is None:
        base.update({"status": "HOLD_COMPONENT_NOT_AVAILABLE" if spec["kind"] in {"MAPPED", "MIXED"} else spec["authority"], "hold_reason": f"shared chain {spec['chain']} output is missing"})
        return base
    if chain_status != "COMPUTED":
        if spec["kind"] == "MIXED":
            aggregate = _aggregate_projection(chain)
            validation = _expected_validation(object_id, aggregate)
            rows = aggregate.get("rows", []) if isinstance(aggregate, Mapping) else []
            primary_present = any(
                str(row.get("cohort")) == "GSE78220" and str(row.get("status")) == "COMPUTED"
                for row in rows if isinstance(row, Mapping)
            )
            if primary_present and validation.get("status") == "PASS":
                base.update({"status": spec["authority"], "aggregate": aggregate, "expected_validation": validation, "hold_reason": "primary D cohort computed; comparison-only cohorts remain non-primary"})
                return base
            if primary_present and validation.get("status") == "MISMATCH":
                base.update({"status": "HOLD_EXPECTED_VALUE_MISMATCH", "aggregate": aggregate, "expected_validation": validation, "hold_reason": "canonical aggregate did not match the frozen output contract"})
                return base
        if spec["chain"] == "D" and spec["kind"] == "HOLD":
            aggregate = _aggregate_projection(chain)
            validation = _expected_validation(object_id, aggregate)
            base.update({"status": "HOLD_EXPECTED_VALUE_MISMATCH" if validation.get("status") == "MISMATCH" else "HOLD_PARTIAL_COHORTS", "aggregate": aggregate, "expected_validation": validation if validation.get("checks") else {"status": "NOT_RUN", "reason": "partial cohort aggregate unavailable"}, "hold_reason": chain.get("hold_reason", "declared D cohort set is partial")})
            return base
        base.update({"status": "HOLD_COMPONENT_NOT_AVAILABLE" if spec["kind"] == "MAPPED" else ("HOLD_PARTIAL_COHORTS" if spec["chain"] == "D" else spec["authority"]), "hold_reason": chain.get("hold_reason", chain.get("reason", f"shared chain {spec['chain']} is held"))})
        return base
    aggregate = _b_projection(chain) if spec["chain"] == "B" else _c_projection(chain) if spec["chain"] == "C" else _aggregate_projection(chain)
    validation = _expected_validation(object_id, aggregate)
    if spec["kind"] == "MIXED" and validation["status"] == "PASS":
        base.update({"status": spec["authority"], "aggregate": aggregate, "expected_validation": validation, "qualification_note": "primary D cohort computed; comparison-only cohorts remain non-primary"})
    elif spec["kind"] == "MIXED" and validation["status"] == "MISMATCH":
        base.update({"status": "HOLD_EXPECTED_VALUE_MISMATCH", "aggregate": aggregate, "expected_validation": validation, "hold_reason": "canonical aggregate did not match the frozen output contract"})
    elif spec["kind"] == "HOLD":
        validation = _expected_validation(object_id, aggregate)
        base.update({"status": "HOLD_PARTIAL_COHORTS" if spec["chain"] == "D" else spec["authority"], "aggregate": aggregate, "expected_validation": validation if validation["checks"] else {"status": "NOT_RUN", "reason": "declared component is not available"}, "hold_reason": "declared adapter component is not available"})
    elif validation["status"] != "PASS":
        base.update({"status": "HOLD_EXPECTED_VALUE_MISMATCH", "aggregate": aggregate, "expected_validation": validation, "hold_reason": "canonical aggregate did not match the frozen output contract"})
    else:
        base.update({"status": "COMPUTED", "aggregate": aggregate, "expected_validation": validation})
    return base


def run_object_adapter(object_id: str, *, shared_dir: str | Path, output_dir: str | Path, repo_root: str | Path | None = None) -> int:
    root = Path(repo_root) if repo_root is not None else Path(__file__).resolve().parents[2]
    payload = build_adapter_payload(object_id, shared_dir=shared_dir, repo_root=root)
    destination = Path(output_dir) / f"{object_id}.status.json"
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n", encoding="utf-8")
    print(f"{object_id}: {payload['status']}")
    return 0 if payload["status"] in {"COMPUTED", "COMPUTED_PRIMARY_PLUS_COMPARISON_ONLY", "CONCEPTUAL_ASSET_CONTRACT", "EXCLUDED_BY_POLICY", "RETIRED"} else 2


def all_object_ids(repo_root: str | Path | None = None) -> list[str]:
    root = repo_root or Path(__file__).resolve().parents[2]
    return list(load_object_map(root))
