"""Frozen IFNG6 paired information decomposition audit.

The module keeps the fail-closed input gate, preregistered decision contract,
and output schema together because they form one lineage-sensitive audit path.
Residual and held-out model fitting are separated into
``ifng6_decomposition_models``; further splitting would detach the decision
fields from the exact machine outputs they adjudicate.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

from .falsification_first import _git_value, _sha256, _write_csv, _write_json
from .falsification_statistics import clustered_auc_analysis, clustered_correlation_analysis
from .ifng6_decomposition_models import lopo_incremental_information, residualized_change


ANALYSIS_ID = "IFNG6_PAIRED_INFORMATION_DECOMPOSITION_V1"
BOOTSTRAP_REPEATS = 10000
PRIMARY_SEED = 20260808
INFORMATIVE_AUC_MIN = 0.70
MATERIAL_AUC_GAP = 0.10
SIGN_DIRECTION_SOURCE = (
    "project-locked IFNG6 equal-gene mean; higher expression is the prespecified "
    "immune-response orientation inherited from the prior audit; no label-informed sign inversion"
)


def _scope_contract() -> str:
    return """# IFNG6 paired information decomposition: scope and frozen contract

This audit directly consumes the prior falsification-first v2 machine outputs. It
does not recompute IFNG6 expression, change the 15-patient/16-record population,
alter response labels, or promote a sensitivity result into a canonical result.

This corrected contract supersedes a retained preliminary run whose decision rule
used therapy-record permutation P despite the required patient-cluster resampling
unit. The correction follows the original unit-of-analysis requirement; it is
recorded as a post-preliminary amendment rather than represented as preregistration.

## Fixed estimands

- PRE, ON, and RAW_DELTA use the prior IFNG6 score table without sign inversion.
- RAW_DELTA = ON - PRE.
- Full residualization fits ON ~ 1 + PRE without response labels.
- Primary residualization is leave-one-patient-out; all records from a held-out
  patient are excluded from its residual fit.
- AUC direction is fixed: higher IFNG6 or positive residual is the positive score.
- Primary uncertainty uses 10,000 deterministic patient-cluster bootstrap draws
  with seed 20260808. Patient 13 remains one resampling cluster.
- Prior exact therapy-record permutation P values are inherited only as
  sensitivities. A patient-cluster label permutation is not defined because
  patient 13 has regimen-specific discordant labels, so no new permutation is run.
- AUCs are calculated directly. No covariance identity is treated as an AUC identity.

## Frozen interpretation rules

An informative channel requires direction-prespecified AUC >= 0.70 and a
patient-cluster bootstrap 95% CI lower bound > 0.50. A material AUC gap is a
point difference >= 0.10; its paired patient-cluster interval is always reported
and controls the strength of wording. Prior record-label exact P values are
inherited as sensitivity results for PRE, ON, and RAW_DELTA only; they do not
adjudicate this decomposition. No residual-channel label permutation is run.

Decision precedence:

1. RAW_CHANGE_ARTIFACT_OR_PARAMETERIZATION_LIMITATION: RAW_DELTA is not
   informative but LOPO residualized change is informative.
2. STATE_MOVEMENT_DISSOCIATION_NOT_SUPPORTED_FOR_IFNG6: PRE, ON, and RAW_DELTA
   are all informative.
3. STATE_INFORMATION_PERSISTS_CHANGE_NOT_READOUT: PRE and ON are informative;
   RAW_DELTA and LOPO residualized change are not; both state channels exceed both
   change channels by at least 0.10 AUC.
4. ON_TREATMENT_INFORMATION_NOT_CAPTURED_BY_DISPLACEMENT: PRE is not informative,
   ON is informative, both change channels are not, and ON exceeds both by at
   least 0.10 AUC.
5. DECOMPOSITION_UNRESOLVED: none of the prespecified patterns is satisfied.

The incremental analysis uses fixed low-dimensional L2 logistic models under
patient-LOPO. It is declared informative only if both the held-out AUC-difference
and equal-patient-weight log-loss improvement intervals exclude zero in the
favorable direction; otherwise it is INCREMENTAL_ANALYSIS_UNSTABLE.
"""


def _read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def _verify_manifest(audit_dir: Path) -> dict[str, str]:
    manifest_path = audit_dir / "11_OUTPUT_MANIFEST.sha256"
    declared = {}
    for line in manifest_path.read_text(encoding="ascii").splitlines():
        digest, name = line.split("  ", 1)
        declared[name] = digest
    for name, digest in declared.items():
        path = audit_dir / name
        if not path.is_file() or _sha256(path) != digest:
            raise ValueError(f"prior audit manifest mismatch: {name}")
    return declared


def _load_frozen_input(prior_audit_dir: Path) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    manifest = _verify_manifest(prior_audit_dir)
    prior_final = _read_json(prior_audit_dir / "10_FINAL_DECISION.json")
    overlap = _read_json(prior_audit_dir / "02_DERIVATION_LONGITUDINAL_OVERLAP.json")
    provenance = _read_json(prior_audit_dir / "01_PROVENANCE.json")
    if prior_final.get("CENTRAL_DECISION") != "GENERAL_STORY_SUPPORTED_WITH_NARROWING":
        raise ValueError("prior central decision does not match the inherited audit state")
    if int(overlap.get("N_overlap_patients", -1)) != 0:
        raise ValueError("derivation/longitudinal overlap is not zero")

    score_path = prior_audit_dir / "04_PAIRED_SIGNATURE_SCORES.csv"
    with score_path.open(encoding="utf-8", newline="") as handle:
        source_rows = [row for row in csv.DictReader(handle) if row["signature"] == "IFNG6"]
    rows = []
    for row in source_rows:
        pre = float(row["PRE_score"])
        on = float(row["ON_score"])
        raw_delta = float(row["delta_score"])
        if not math.isclose(raw_delta, on - pre, rel_tol=0.0, abs_tol=1e-12):
            raise ValueError(f"prior IFNG6 delta identity mismatch: {row['record_id']}")
        rows.append({
            "patient_id": str(row["patient_id"]),
            "therapy_record_id": row["record_id"],
            "accession": "PRJEB23709",
            "response": row["response_label"],
            "response_binary": int(row["response_binary"]),
            "regimen": row["therapy_type"],
            "IFNG6_PRE": pre,
            "IFNG6_ON": on,
            "IFNG6_RAW_DELTA": raw_delta,
        })
    patients = set(row["patient_id"] for row in rows)
    record_ids = [row["therapy_record_id"] for row in rows]
    if len(rows) != 16 or len(patients) != 15 or len(record_ids) != len(set(record_ids)):
        raise ValueError(f"frozen paired population conflict: {len(rows)} records/{len(patients)} patients")
    patient_13 = [row for row in rows if row["patient_id"] == "13"]
    if len(patient_13) != 2:
        raise ValueError("patient 13 must retain both therapy records")

    matrix_path = prior_audit_dir / "05_STATE_MOVEMENT_MATRIX.csv"
    with matrix_path.open(encoding="utf-8", newline="") as handle:
        locked_rows = [row for row in csv.DictReader(handle) if row["signature"] == "IFNG6"]
    if len(locked_rows) != 1:
        raise ValueError("prior IFNG6 matrix row is missing or duplicated")
    return rows, {
        "manifest": manifest,
        "prior_final": prior_final,
        "overlap": overlap,
        "provenance": provenance,
        "locked_ifng6": locked_rows[0],
    }


def _assert_locked_regression(primary: Mapping[str, Any], locked: Mapping[str, str]) -> None:
    mapping = {
        "PRE": ("PRE", "PRE_AUC", "PRE_CI_low", "PRE_CI_high"),
        "ON": ("ON", "ON_AUC", "ON_CI_low", "ON_CI_high"),
        "RAW_DELTA": (
            "RAW_DELTA", "MOVEMENT_AUC", "MOVEMENT_CI_low", "MOVEMENT_CI_high"
        ),
    }
    for _, (channel, auc_key, low_key, high_key) in mapping.items():
        observed = primary["channels"][channel]
        checks = (
            (observed["auc"], float(locked[auc_key]), auc_key),
            (observed["auc_ci_95"][0], float(locked[low_key]), low_key),
            (observed["auc_ci_95"][1], float(locked[high_key]), high_key),
        )
        for current, expected, name in checks:
            if not math.isclose(float(current), expected, rel_tol=0.0, abs_tol=1e-12):
                raise ValueError(f"prior IFNG6 regression mismatch for {name}: {current} != {expected}")
    comparison = primary["comparisons"]["ON_MINUS_RAW_DELTA"]
    checks = (
        (comparison["estimate"], float(locked["on_minus_movement_auc"]), "on_minus_movement_auc"),
        (comparison["ci_95"][0], float(locked["on_minus_movement_ci_low"]), "on_minus_movement_ci_low"),
        (comparison["ci_95"][1], float(locked["on_minus_movement_ci_high"]), "on_minus_movement_ci_high"),
    )
    for current, expected, name in checks:
        if not math.isclose(float(current), expected, rel_tol=0.0, abs_tol=1e-12):
            raise ValueError(f"prior IFNG6 regression mismatch for {name}: {current} != {expected}")


def _is_informative(channel: Mapping[str, Any]) -> bool:
    return float(channel["auc"]) >= INFORMATIVE_AUC_MIN and float(channel["auc_ci_95"][0]) > 0.5


def _material_gap(results: Mapping[str, Any], left: str, right: str) -> bool:
    return float(results["comparisons"][f"{left}_MINUS_{right}"]["estimate"]) >= MATERIAL_AUC_GAP


def _decomposition_decision(results: Mapping[str, Any]) -> str:
    channels = results["channels"]
    pre = _is_informative(channels["PRE"])
    on = _is_informative(channels["ON"])
    raw = _is_informative(channels["RAW_DELTA"])
    residual = _is_informative(channels["RESIDUAL_LOPO"])
    if not raw and residual:
        return "RAW_CHANGE_ARTIFACT_OR_PARAMETERIZATION_LIMITATION"
    if pre and on and raw:
        return "STATE_MOVEMENT_DISSOCIATION_NOT_SUPPORTED_FOR_IFNG6"
    if pre and on and not raw and not residual and all(
        _material_gap(results, state, change)
        for state in ("PRE", "ON") for change in ("RAW_DELTA", "RESIDUAL_LOPO")
    ):
        return "STATE_INFORMATION_PERSISTS_CHANGE_NOT_READOUT"
    if not pre and on and not raw and not residual and all(
        _material_gap(results, "ON", change) for change in ("RAW_DELTA", "RESIDUAL_LOPO")
    ):
        return "ON_TREATMENT_INFORMATION_NOT_CAPTURED_BY_DISPLACEMENT"
    return "DECOMPOSITION_UNRESOLVED"


def _propositions(decision: str) -> tuple[str, list[str]]:
    allowed = {
        "STATE_INFORMATION_PERSISTS_CHANGE_NOT_READOUT": "Response information is predominantly state-like and already present before or independent of the measured displacement; treatment-associated IFNG6 change should not automatically be interpreted as a benefit-relevant pharmacodynamic readout.",
        "ON_TREATMENT_INFORMATION_NOT_CAPTURED_BY_DISPLACEMENT": "Response information is enriched in the on-treatment state but is not represented by simple longitudinal displacement; nonlinear transition, timing, composition, saturation, and measurement geometry remain alternatives.",
        "RAW_CHANGE_ARTIFACT_OR_PARAMETERIZATION_LIMITATION": "Raw difference-score parameterization obscured response-relevant baseline-adjusted change; raw delta null evidence cannot support the central state-movement claim for IFNG6.",
        "STATE_MOVEMENT_DISSOCIATION_NOT_SUPPORTED_FOR_IFNG6": "State-movement dissociation is not supported for IFNG6 because state and change channels all carry response information.",
        "DECOMPOSITION_UNRESOLVED": "The current paired cohort does not uniquely allocate IFNG6 response information between state and baseline-adjusted change.",
    }[decision]
    forbidden = [
        "treatment has no effect",
        "movement predicts the opposite response",
        "the residualized score identifies a causal treatment effect",
        "AUC obeys an additive PRE plus change identity",
        "the 15-patient decomposition establishes a universal dynamics law",
    ]
    return allowed, forbidden


def _final_report(decision: Mapping[str, Any], results: Mapping[str, Any], incremental: Mapping[str, Any]) -> str:
    channels = results["channels"]
    lines = [
        "# IFNG6 paired information decomposition",
        "",
        f"FINAL_DECOMPOSITION_DECISION = {decision['FINAL_DECOMPOSITION_DECISION']}",
        "",
        "| Channel | Direction-fixed AUC | Patient-cluster 95% CI | Inherited record-label P |",
        "|---|---:|---:|---:|",
    ]
    for name in ("PRE", "ON", "RAW_DELTA", "RESIDUAL_FULL", "RESIDUAL_LOPO"):
        channel = channels[name]
        ci = channel["auc_ci_95"]
        inherited_p = decision["INHERITED_RECORD_LABEL_P_SENSITIVITY"].get(name)
        p_text = "NA" if inherited_p is None else f"{inherited_p:.4f}"
        lines.append(f"| {name} | {channel['auc']:.3f} | [{ci[0]:.3f}, {ci[1]:.3f}] | {p_text} |")
    lines.extend([
        "",
        f"Incremental ON given PRE: {incremental['incremental_on_given_pre']['status']}.",
        f"Incremental RAW_DELTA given ON: {incremental['incremental_delta_given_on']['status']}.",
        "",
        decision["CENTRAL_PROPOSITION_ALLOWED"],
        "",
        "Intervals are patient-cluster intervals over 15 patients and 16 therapy records. They do not separate measurement error, biopsy composition, or regression to the mean.",
    ])
    return "\n".join(lines) + "\n"


def run_ifng6_information_decomposition(
    *,
    prior_audit_dir: str | Path,
    repository_root: str | Path,
    output_dir: str | Path,
    bootstrap_repeats: int = BOOTSTRAP_REPEATS,
    seed: int = PRIMARY_SEED,
) -> dict[str, Any]:
    prior = Path(prior_audit_dir).resolve()
    repository = Path(repository_root).resolve()
    output = Path(output_dir).resolve()
    if bootstrap_repeats != BOOTSTRAP_REPEATS or seed != PRIMARY_SEED:
        raise ValueError(f"frozen decomposition requires {BOOTSTRAP_REPEATS} bootstrap draws and seed {PRIMARY_SEED}")
    if output.exists() and any(output.iterdir()):
        raise ValueError(f"output directory must be absent or empty: {output}")
    output.mkdir(parents=True, exist_ok=True)
    (output / "00_SCOPE_AND_CONTRACT.md").write_text(_scope_contract(), encoding="utf-8")

    rows, inherited = _load_frozen_input(prior)
    pre = [float(row["IFNG6_PRE"]) for row in rows]
    on = [float(row["IFNG6_ON"]) for row in rows]
    raw_delta = [float(row["IFNG6_RAW_DELTA"]) for row in rows]
    labels = [int(row["response_binary"]) for row in rows]
    patient_ids = [str(row["patient_id"]) for row in rows]

    residual = residualized_change(pre, on, patient_ids)
    residual_full = residual["full_sample"]["residuals"]
    residual_lopo = residual["lopo"]["residuals"]
    channels = {
        "PRE": pre,
        "ON": on,
        "RAW_DELTA": raw_delta,
        "RESIDUAL_FULL": residual_full,
        "RESIDUAL_LOPO": residual_lopo,
    }
    comparisons = [
        ("ON", "PRE"),
        ("ON", "RAW_DELTA"),
        ("PRE", "RAW_DELTA"),
        ("ON", "RESIDUAL_FULL"),
        ("PRE", "RESIDUAL_FULL"),
        ("RAW_DELTA", "RESIDUAL_FULL"),
        ("ON", "RESIDUAL_LOPO"),
        ("PRE", "RESIDUAL_LOPO"),
        ("RAW_DELTA", "RESIDUAL_LOPO"),
    ]
    auc_results, bootstrap_rows = clustered_auc_analysis(
        channels, labels, patient_ids, comparisons,
        bootstrap_repeats=bootstrap_repeats, seed=seed,
        compute_record_label_permutation=False,
    )
    _assert_locked_regression(auc_results, inherited["locked_ifng6"])
    for row, full_value, lopo_value in zip(rows, residual_full, residual_lopo):
        row["IFNG6_RESIDUAL_FULL"] = full_value
        row["IFNG6_RESIDUAL_LOPO"] = lopo_value

    mean_state = ((np.asarray(pre) + np.asarray(on)) / 2.0).tolist()
    diagnostics = {
        "PRE_vs_RAW_DELTA": clustered_correlation_analysis(
            pre, raw_delta, patient_ids, bootstrap_repeats=bootstrap_repeats, seed=seed + 100
        ),
        "MEAN_STATE_vs_RAW_DELTA": clustered_correlation_analysis(
            mean_state, raw_delta, patient_ids, bootstrap_repeats=bootstrap_repeats, seed=seed + 101
        ),
        "interpretation_boundary": "PRE-delta association is mathematically coupled; without technical replicates, measurement error and regression to the mean cannot be separated from biology.",
    }
    incremental = lopo_incremental_information(rows, bootstrap_repeats=bootstrap_repeats, seed=seed + 200)
    decomposition = _decomposition_decision(auc_results)
    allowed, forbidden = _propositions(decomposition)
    comparisons_out = auc_results["comparisons"]
    final = {
        "N_PATIENTS": auc_results["n_patients"],
        "N_RECORDS": auc_results["n_records"],
        "PRE_AUC": auc_results["channels"]["PRE"]["auc"],
        "PRE_AUC_CI": auc_results["channels"]["PRE"]["auc_ci_95"],
        "ON_AUC": auc_results["channels"]["ON"]["auc"],
        "ON_AUC_CI": auc_results["channels"]["ON"]["auc_ci_95"],
        "RAW_DELTA_AUC": auc_results["channels"]["RAW_DELTA"]["auc"],
        "RAW_DELTA_AUC_CI": auc_results["channels"]["RAW_DELTA"]["auc_ci_95"],
        "RESIDUAL_CHANGE_AUC_FULL": auc_results["channels"]["RESIDUAL_FULL"]["auc"],
        "RESIDUAL_CHANGE_AUC_LOPO": auc_results["channels"]["RESIDUAL_LOPO"]["auc"],
        "RESIDUAL_CHANGE_AUC_CI": auc_results["channels"]["RESIDUAL_LOPO"]["auc_ci_95"],
        "ON_MINUS_PRE_AUC": comparisons_out["ON_MINUS_PRE"]["estimate"],
        "ON_MINUS_PRE_CI": comparisons_out["ON_MINUS_PRE"]["ci_95"],
        "ON_MINUS_RAW_DELTA_AUC": comparisons_out["ON_MINUS_RAW_DELTA"]["estimate"],
        "ON_MINUS_RAW_DELTA_CI": comparisons_out["ON_MINUS_RAW_DELTA"]["ci_95"],
        "ON_MINUS_RESIDUAL_CHANGE_AUC": comparisons_out["ON_MINUS_RESIDUAL_LOPO"]["estimate"],
        "ON_MINUS_RESIDUAL_CHANGE_CI": comparisons_out["ON_MINUS_RESIDUAL_LOPO"]["ci_95"],
        "COR_PRE_RAW_DELTA": diagnostics["PRE_vs_RAW_DELTA"]["spearman_rho"],
        "COR_MEAN_RAW_DELTA": diagnostics["MEAN_STATE_vs_RAW_DELTA"]["spearman_rho"],
        "INCREMENTAL_ON_GIVEN_PRE": incremental["incremental_on_given_pre"]["status"],
        "INCREMENTAL_DELTA_GIVEN_ON": incremental["incremental_delta_given_on"]["status"],
        "SIGN_DIRECTION_SOURCE": SIGN_DIRECTION_SOURCE,
        "INHERITED_RECORD_LABEL_P_SENSITIVITY": {
            "PRE": float(inherited["locked_ifng6"]["PRE_permutation_p"]),
            "ON": float(inherited["locked_ifng6"]["ON_permutation_p"]),
            "RAW_DELTA": float(inherited["locked_ifng6"]["MOVEMENT_permutation_p"]),
            "RESIDUAL_FULL": None,
            "RESIDUAL_LOPO": None,
        },
        "FINAL_DECOMPOSITION_DECISION": decomposition,
        "CENTRAL_PROPOSITION_ALLOWED": allowed,
        "CENTRAL_PROPOSITION_FORBIDDEN": forbidden,
        "STRONGEST_UNRESOLVED_ALTERNATIVE": "baseline-adjusted change, measurement error, biopsy composition, nonlinear timing, and regression to the mean remain incompletely separated in one 15-patient cohort",
        "SUPERSEDES_PRELIMINARY_REASON": "the preliminary decomposition incorrectly used therapy-record label permutation P as a decision threshold despite the required patient-cluster resampling unit",
    }

    primary_names = ("PRE", "ON", "RAW_DELTA")
    primary_results = {
        **{key: auc_results[key] for key in ("n_records", "n_patients", "n_R", "n_NR", "bootstrap_repeats_requested", "bootstrap_repeats_valid", "permutation_unit_note")},
        "channels": {name: auc_results["channels"][name] for name in primary_names},
        "comparisons": {
            key: value for key, value in comparisons_out.items()
            if key in {"ON_MINUS_PRE", "ON_MINUS_RAW_DELTA", "PRE_MINUS_RAW_DELTA"}
        },
        "prior_locked_regression": "PASS",
        "inherited_record_label_permutation_p_sensitivity": {
            "PRE": float(inherited["locked_ifng6"]["PRE_permutation_p"]),
            "ON": float(inherited["locked_ifng6"]["ON_permutation_p"]),
            "RAW_DELTA": float(inherited["locked_ifng6"]["MOVEMENT_permutation_p"]),
        },
    }
    residual_results = {
        "model": residual,
        "channels": {name: auc_results["channels"][name] for name in ("RESIDUAL_FULL", "RESIDUAL_LOPO")},
        "comparisons": {key: value for key, value in comparisons_out.items() if "RESIDUAL" in key},
        "primary_version": "RESIDUAL_LOPO",
        "response_labels_used_in_residualization": False,
    }
    provenance = {
        "analysis_id": ANALYSIS_ID,
        "repository_root": str(repository),
        "repository_head": _git_value(repository, "rev-parse", "HEAD"),
        "repository_dirty": bool(_git_value(repository, "status", "--porcelain")),
        "prior_audit_dir": str(prior),
        "prior_analysis_id": inherited["provenance"]["analysis_id"],
        "prior_central_decision": inherited["prior_final"]["CENTRAL_DECISION"],
        "derivation_longitudinal_overlap": inherited["overlap"]["N_overlap_patients"],
        "parameters": {"bootstrap_repeats": bootstrap_repeats, "primary_seed": seed, "incremental_seed": seed + 200},
        "decision_thresholds": {"informative_auc_min": INFORMATIVE_AUC_MIN, "patient_cluster_ci_lower_min_exclusive": 0.5, "material_auc_gap": MATERIAL_AUC_GAP},
        "contract_amendment": {
            "reason": "patient 13 has discordant regimen-specific labels, so therapy-record label permutation is not a valid patient-cluster primary inference",
            "timing": "amended after preliminary run exposed the resampling-unit conflict; preliminary output retained and not promoted",
            "preliminary_output_dir": str(output.parent / output.name.removesuffix("_v2")),
        },
        "inputs": {
            name: {"sha256": digest} for name, digest in inherited["manifest"].items()
            if name in {"01_PROVENANCE.json", "02_DERIVATION_LONGITUDINAL_OVERLAP.json", "04_PAIRED_SIGNATURE_SCORES.csv", "05_STATE_MOVEMENT_MATRIX.csv", "10_FINAL_DECISION.json"}
        },
        "implementation": {
            name: {"path": path, "sha256": _sha256(repository / path)}
            for name, path in {
                "runner": "scripts/run_ifng6_information_decomposition.py",
                "orchestration": "src/meld_icb/ifng6_information_decomposition.py",
                "models": "src/meld_icb/ifng6_decomposition_models.py",
                "statistics": "src/meld_icb/falsification_statistics.py",
            }.items()
        },
        "sign_direction_source": SIGN_DIRECTION_SOURCE,
    }

    _write_json(output / "01_INPUT_PROVENANCE.json", provenance)
    _write_csv(output / "02_FROZEN_PAIRED_DATA.csv", rows)
    _write_json(output / "03_PRE_ON_DELTA_RESULTS.json", primary_results)
    _write_json(output / "04_RESIDUALIZED_CHANGE_RESULTS.json", residual_results)
    _write_json(output / "05_INCREMENTAL_INFORMATION.json", incremental)
    _write_json(output / "06_CHANGE_SCORE_DIAGNOSTICS.json", diagnostics)
    _write_csv(output / "07_PAIRED_BOOTSTRAP_RESULTS.csv", bootstrap_rows)
    _write_json(output / "08_INTERPRETATION_DECISION.json", final)
    (output / "09_FINAL_REPORT.md").write_text(_final_report(final, auc_results, incremental), encoding="utf-8")
    manifest_rows = [
        f"{_sha256(path)}  {path.name}"
        for path in sorted(output.iterdir(), key=lambda item: item.name)
        if path.is_file() and path.name != "10_OUTPUT_MANIFEST.sha256"
    ]
    (output / "10_OUTPUT_MANIFEST.sha256").write_text("\n".join(manifest_rows) + "\n", encoding="ascii")
    return {
        "status": "COMPUTED",
        "output_dir": str(output),
        "final_decomposition_decision": decomposition,
    }
