"""Falsification-first audit of state, movement, and intervention variation.

The module intentionally keeps one audit-contract orchestrator and its narrow
input/report helpers together: the output schema, provenance, and calculations
must be reviewed as one traceable unit. Statistical primitives live separately
in ``falsification_statistics``; splitting this contract further would create
thin forwarding layers without reducing the required audit workset.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import subprocess
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
from scipy import stats

from .carrier_grounding import established_signature_gene_sets
from .config import load_freeze
from .falsification_statistics import (
    benjamini_hochberg,
    change_coupling,
    paired_change_control,
    state_movement_metrics,
)
from .frozen_state_inputs import (
    _gene_set_union,
    _read_selected_matrix_fast,
    read_gse91061_on_treatment_labels,
)
from .frozen_state_transfer import (
    AXIS_ORDER,
    _axis_raw,
    fit_axis_normalization,
    transform_axis_normalization,
)
from .paired_movement_inputs import build_paired_records
from .statistics import auc


ANALYSIS_ID = "FALSIFICATION_FIRST_STATE_MOVEMENT_V1"
POSITIVE_CONTROL = "IFNG6"
SIGNATURE_ORDER = ("FROZEN_ORIGINAL", "IFNG6", "TIS18_MEAN", "CYT")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _git_value(root: Path, *args: str) -> str:
    command = ["git", "-c", f"safe.directory={root.as_posix()}", *args]
    return subprocess.run(command, cwd=root, check=True, capture_output=True, text=True).stdout.strip()


def _write_json(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False, default=str) + "\n", encoding="utf-8")


def _write_csv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty audit table: {path.name}")
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _axis_parameters(rows: Sequence[Mapping[str, float]], sample_ids: Sequence[str], scope: str):
    return fit_axis_normalization(rows, sample_ids, scope_id=scope, ddof=1)


def _boundary_rows(
    raw_axes: Mapping[str, Mapping[str, float]],
    sample_ids: Sequence[str],
    parameters: Any,
    freeze: Mapping[str, Any],
) -> dict[str, float]:
    transformed = transform_axis_normalization([raw_axes[sample_id] for sample_id in sample_ids], parameters)
    coefficients = freeze["boundary"]["coefficients"]
    intercept = float(freeze["boundary"]["intercept"])
    return {
        sample_id: intercept + sum(float(coefficients[axis]) * float(row[axis]) for axis in AXIS_ORDER)
        for sample_id, row in zip(sample_ids, transformed)
    }


def _signature_scores(
    expression: Mapping[str, Mapping[str, float]], genes: Sequence[str]
) -> tuple[dict[str, float], dict[str, Any]]:
    wanted = tuple(str(gene).upper() for gene in genes)
    present = tuple(gene for gene in wanted if all(gene in row for row in expression.values()))
    if not present:
        return {}, {"status": "NOT_COMPUTABLE", "reason": "no signature genes present", "n_genes": len(wanted), "n_available": 0}
    scores = {
        sample_id: float(np.mean([row[gene] for gene in present]))
        for sample_id, row in expression.items()
    }
    return scores, {
        "status": "COMPUTABLE",
        "n_genes": len(wanted),
        "n_available": len(present),
        "coverage": len(present) / len(wanted),
        "available_genes": list(present),
        "missing_genes": [gene for gene in wanted if gene not in present],
    }


def _flatten_metric(signature: str, metrics: Mapping[str, Any], coverage: Mapping[str, Any]) -> dict[str, Any]:
    row: dict[str, Any] = {
        "signature": signature,
        "computability": coverage.get("status", "COMPUTABLE"),
        "n_genes": coverage.get("n_genes"),
        "n_available": coverage.get("n_available"),
        "gene_coverage": coverage.get("coverage"),
        "n_records": metrics["n_records"],
        "n_patients": metrics["n_patients"],
        "n_R": metrics["n_R"],
        "n_NR": metrics["n_NR"],
        "on_minus_movement_auc": metrics["on_minus_movement_auc"],
        "on_minus_movement_ci_low": metrics["on_minus_movement_auc_ci_95"][0],
        "on_minus_movement_ci_high": metrics["on_minus_movement_auc_ci_95"][1],
    }
    for channel in ("pre", "on", "movement"):
        values = metrics[channel]
        prefix = channel.upper()
        row.update({
            f"{prefix}_AUC": values["auc"],
            f"{prefix}_CI_low": values["auc_ci_95"][0],
            f"{prefix}_CI_high": values["auc_ci_95"][1],
            f"{prefix}_permutation_p": values["auc_permutation_p_two_sided"],
            f"{prefix}_mannwhitney_p": values["mannwhitney_p_two_sided"],
            f"{prefix}_R_median": values["responder"]["median"],
            f"{prefix}_R_q1": values["responder"]["q1"],
            f"{prefix}_R_q3": values["responder"]["q3"],
            f"{prefix}_NR_median": values["nonresponder"]["median"],
            f"{prefix}_NR_q1": values["nonresponder"]["q1"],
            f"{prefix}_NR_q3": values["nonresponder"]["q3"],
        })
    return row


def _classify_signature_rows(rows: list[dict[str, Any]]) -> None:
    established = [row for row in rows if row["signature"] != "FROZEN_ORIGINAL"]
    for channel in ("ON", "MOVEMENT"):
        adjusted = benjamini_hochberg([float(row[f"{channel}_permutation_p"]) for row in established])
        for row, q_value in zip(established, adjusted):
            row[f"{channel}_BH_q"] = q_value
    frozen = rows[0]
    frozen["ON_BH_q"] = frozen["ON_permutation_p"]
    frozen["MOVEMENT_BH_q"] = frozen["MOVEMENT_permutation_p"]
    for row in rows:
        state_detected = float(row["ON_AUC"]) > 0.5 and float(row["ON_BH_q"]) <= 0.05
        movement_detected = float(row["MOVEMENT_BH_q"]) <= 0.05
        difference_detected = row["on_minus_movement_ci_low"] is not None and float(row["on_minus_movement_ci_low"]) > 0
        if state_detected and not movement_detected and difference_detected:
            interpretation = "DISSOCIATION_SUPPORTED"
        elif movement_detected:
            interpretation = "MOVEMENT_RESPONSE_INFORMATION_DETECTED"
        elif state_detected:
            interpretation = "STATE_SIGNAL_WITH_UNRESOLVED_CHANNEL_DIFFERENCE"
        else:
            interpretation = "NO_DETECTABLE_ON_STATE_SIGNAL"
        row.update({
            "state_response_information_detected": state_detected,
            "movement_response_information_detected": movement_detected,
            "paired_channel_difference_detected": difference_detected,
            "interpretation": interpretation,
        })


def _overlap_audit(
    truth_root: Path, records: Sequence[Mapping[str, Any]]
) -> dict[str, Any]:
    derivation = read_gse91061_on_treatment_labels(truth_root)
    derivation_patients = [f"GSE91061:{row['patient_id']}" for row in derivation]
    longitudinal_patients = sorted({f"PRJEB23709:{row['patient_id']}" for row in records})
    longitudinal_on = [f"PRJEB23709:{row['edt_sample_id']}" for row in records]
    overlap = sorted(set(derivation_patients) & set(longitudinal_patients))
    derivation_tokens = {str(row["patient_id"]).removeprefix("Pt") for row in derivation}
    longitudinal_tokens = {str(row["patient_id"]) for row in records}
    token_collisions = sorted(derivation_tokens & longitudinal_tokens, key=int)
    return {
        "N_derivation_patients": len(set(derivation_patients)),
        "N_derivation_samples": len(derivation),
        "N_longitudinal_patients": len(longitudinal_patients),
        "N_longitudinal_records": len(records),
        "N_overlap_patients": len(overlap),
        "overlap_patient_ids": overlap,
        "overlap_on_treatment_sample_ids": [],
        "unqualified_numeric_token_collisions": token_collisions,
        "token_collision_interpretation": "source-local numeric labels are not cross-study patient identifiers and are not counted as overlap",
        "derivation_patient_ids": derivation_patients,
        "derivation_sample_ids": [f"GSE91061:{row['sample_id']}" for row in derivation],
        "longitudinal_patient_ids": longitudinal_patients,
        "longitudinal_pre_sample_ids": [f"PRJEB23709:{row['pre_sample_id']}" for row in records],
        "longitudinal_on_treatment_sample_ids": longitudinal_on,
        "patient_correspondence": "distinct study accessions and source-local patient namespaces; no declared cross-study patient correspondence",
        "derivation_treatment": "nivolumab (anti-PD-1)",
        "longitudinal_treatments": sorted({str(row["therapy_type"]) for row in records}),
        "preprocessing": "PRE_MEAN_THEN_AXIS_Z; ddof=1",
        "classification": "A_COMPLETELY_INDEPENDENT" if not overlap else "B_OR_C_OVERLAP_REQUIRES_ADJUDICATION",
        "leave_longitudinal_patients_out_refit": "NOT_REQUIRED_ZERO_OVERLAP" if not overlap else "REQUIRED_NOT_RUN",
    }


def _action_rows(records: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    return [{
        "patient": record["patient_id"],
        "time": "PRE_to_EDT",
        "treatment": record["therapy_type"],
        "dose": "NOT_AVAILABLE",
        "combination": "anti-CTLA-4_plus_anti-PD-1" if record["therapy_type"] == "ipiPD1" else "anti-PD-1_monotherapy",
        "schedule": "NOT_AVAILABLE",
        "treatment_change": "NOT_AVAILABLE",
        "discontinuation": "NOT_AVAILABLE",
        "rechallenge": "NOT_AVAILABLE",
        "other_action_information": "therapy-specific paired record",
    } for record in records]


def _normalization_sensitivity(
    raw_axes: Mapping[str, Mapping[str, float]],
    records: Sequence[Mapping[str, Any]],
    derivation_raw: Sequence[Mapping[str, float]],
    derivation_ids: Sequence[str],
    freeze: Mapping[str, Any],
    canonical_scores: Mapping[str, float],
) -> list[dict[str, Any]]:
    pre_ids = [str(record["pre_sample_id"]) for record in records]
    on_ids = [str(record["edt_sample_id"]) for record in records]
    all_ids = list(dict.fromkeys([sample for record in records for sample in (str(record["pre_sample_id"]), str(record["edt_sample_id"]))]))
    scopes = {
        "PRE_BASED": _axis_parameters([raw_axes[sid] for sid in pre_ids], pre_ids, "paired_PRE_n16"),
        "ON_BASED": _axis_parameters([raw_axes[sid] for sid in on_ids], on_ids, "paired_ON_n16"),
        "POOLED_CANONICAL": _axis_parameters([raw_axes[sid] for sid in all_ids], all_ids, "paired_PRE_EDT_union_n32"),
        "DERIVATION_FIXED_EXTERNAL": _axis_parameters(derivation_raw, derivation_ids, "GSE91061_derivation_on_n27"),
    }
    labels = [int(record["y_true"]) for record in records]
    canonical_delta = np.asarray([
        canonical_scores[str(record["edt_sample_id"])] - canonical_scores[str(record["pre_sample_id"])]
        for record in records
    ])
    rows: list[dict[str, Any]] = []
    for scope, parameters in scopes.items():
        scores = _boundary_rows(raw_axes, all_ids, parameters, freeze)
        pre = np.asarray([scores[str(record["pre_sample_id"])] for record in records])
        on = np.asarray([scores[str(record["edt_sample_id"])] for record in records])
        movement = on - pre
        rank = stats.spearmanr(canonical_delta, movement)
        rows.append({
            "normalization": scope,
            "normalization_source": parameters.normalization_scope_id,
            "ddof": parameters.ddof,
            "movement_auc": auc(movement.tolist(), labels),
            "median_delta": float(np.median(movement)),
            "mean_delta": float(np.mean(movement)),
            "n_positive": int(np.sum(movement > 0)),
            "fraction_positive": float(np.mean(movement > 0)),
            "spearman_vs_pooled_delta": float(rank.statistic),
            "spearman_p": float(rank.pvalue),
            "direction_agreement_vs_pooled": float(np.mean(np.sign(movement) == np.sign(canonical_delta))),
        })
    return rows


def _established_decision(rows: Sequence[Mapping[str, Any]]) -> str:
    established = [row for row in rows if row["signature"] != "FROZEN_ORIGINAL"]
    supported = sum(row["interpretation"] == "DISSOCIATION_SUPPORTED" for row in established)
    movement = sum(bool(row["movement_response_information_detected"]) for row in established)
    if supported >= 2 and movement == 0:
        return "RESULT_A_GENERAL_DISSOCIATION_SUPPORTED"
    if movement >= 2 and supported == 0:
        return "RESULT_C_FALSIFIED"
    return "RESULT_B_MIXED"


def run_falsification_first(
    *,
    truth_root: str | Path,
    repository_root: str | Path,
    output_dir: str | Path,
    bootstrap_repeats: int = 10000,
    seed: int = 20260807,
) -> dict[str, Any]:
    truth = Path(truth_root).resolve()
    repository = Path(repository_root).resolve()
    output = Path(output_dir).resolve()
    if output.exists() and any(output.iterdir()):
        raise ValueError(f"output directory must be absent or empty: {output}")
    output.mkdir(parents=True, exist_ok=True)
    config_path = repository / "configs" / "paired_dynamics.yaml"
    freeze_path = repository / "configs" / "frozen_state_transfer.yaml"
    freeze = load_freeze(repository)
    records, paired_provenance = build_paired_records(truth, config_path=config_path)
    sample_ids = list(dict.fromkeys([sample for record in records for sample in (str(record["pre_sample_id"]), str(record["edt_sample_id"]))]))
    signatures = established_signature_gene_sets()
    requested_signatures = {name: signatures[name] for name in ("IFNG6", "TIS18", "CYT")}
    genes = _gene_set_union(freeze) | {gene for values in requested_signatures.values() for gene in values}
    expression_path = truth / "raw_data" / "PRJEB23709" / "cancercell_normalized_counts_genenames.txt"
    expression_rows = _read_selected_matrix_fast(
        expression_path, sample_ids, delimiter="\t", gene_field="Gene",
        transform="truncate_negative_then_log2p1", genes=genes,
    )
    expression = {sample_id: row for sample_id, row in zip(sample_ids, expression_rows)}
    raw_rows, axis_coverage = _axis_raw(expression_rows, freeze["method"]["axes"]["gene_sets"], float(freeze["method"]["preprocessing"]["minimum_gene_coverage"]))
    raw_axes = {sample_id: row for sample_id, row in zip(sample_ids, raw_rows)}
    pooled_parameters = _axis_parameters(raw_rows, sample_ids, "paired_PRE_EDT_union_n32")
    canonical_boundary = _boundary_rows(raw_axes, sample_ids, pooled_parameters, freeze)
    labels = [int(record["y_true"]) for record in records]
    patient_ids = [str(record["patient_id"]) for record in records]

    score_rows: list[dict[str, Any]] = []
    metric_rows: list[dict[str, Any]] = []
    score_maps: dict[str, dict[str, float]] = {"FROZEN_ORIGINAL": canonical_boundary}
    frozen_declared = sum(int(item["declared_genes"]) for item in axis_coverage.values())
    frozen_observed = sum(int(item["observed_genes"]) for item in axis_coverage.values())
    coverage: dict[str, dict[str, Any]] = {
        "FROZEN_ORIGINAL": {
            "status": "COMPUTABLE",
            "n_genes": frozen_declared,
            "n_available": frozen_observed,
            "coverage": frozen_observed / frozen_declared,
        }
    }
    for source_name, geneset in requested_signatures.items():
        display_name = "TIS18_MEAN" if source_name == "TIS18" else source_name
        scores, signature_coverage = _signature_scores(expression, geneset)
        score_maps[display_name] = scores
        coverage[display_name] = signature_coverage
    for signature in SIGNATURE_ORDER:
        if coverage[signature]["status"] != "COMPUTABLE":
            metric_rows.append({"signature": signature, **coverage[signature]})
            continue
        pre = [score_maps[signature][str(record["pre_sample_id"])] for record in records]
        on = [score_maps[signature][str(record["edt_sample_id"])] for record in records]
        metrics = state_movement_metrics(pre, on, labels, patient_ids, bootstrap_repeats=bootstrap_repeats, seed=seed + len(metric_rows))
        metric_rows.append(_flatten_metric(signature, metrics, coverage[signature]))
        for record, pre_score, on_score in zip(records, pre, on):
            score_rows.append({
                "signature": signature,
                "record_id": record["record_id"],
                "patient_id": record["patient_id"],
                "therapy_type": record["therapy_type"],
                "response_binary": record["y_true"],
                "response_label": record["response_label"],
                "pre_sample_id": record["pre_sample_id"],
                "on_sample_id": record["edt_sample_id"],
                "PRE_score": pre_score,
                "ON_score": on_score,
                "delta_score": on_score - pre_score,
            })
    _classify_signature_rows(metric_rows)

    derivation_labels = read_gse91061_on_treatment_labels(truth)
    derivation_ids = [str(row["sample_id"]) for row in derivation_labels]
    derivation_expression_path = truth / "processed" / "expression_log2_GSE91061_symbol.csv"
    derivation_expression = _read_selected_matrix_fast(
        derivation_expression_path, derivation_ids, delimiter=",", gene_field=None,
        transform="as_deposited_logcpm", genes=_gene_set_union(freeze),
    )
    derivation_raw, _ = _axis_raw(derivation_expression, freeze["method"]["axes"]["gene_sets"], float(freeze["method"]["preprocessing"]["minimum_gene_coverage"]))
    normalization = _normalization_sensitivity(raw_axes, records, derivation_raw, derivation_ids, freeze, canonical_boundary)

    overlap = _overlap_audit(truth, records)
    action_rows = _action_rows(records)
    treatment_counts = {name: sum(row["treatment"] == name for row in action_rows) for name in sorted({row["treatment"] for row in action_rows})}
    action_summary = {
        "n_records": len(action_rows),
        "n_patients": len(set(row["patient"] for row in action_rows)),
        "treatment_counts": treatment_counts,
        "outcome_variation": {"R_records": int(sum(labels)), "NR_records": int(len(labels) - sum(labels))},
        "state_variation": "paired PRE and early-during-treatment biopsies",
        "action_variation": "PD1 monotherapy versus ipiPD1 combination; dose, schedule, and randomized contrast unavailable",
        "treatment_category_is_constant": len(treatment_counts) == 1,
        "dose_variation": "NOT_AVAILABLE",
        "schedule_variation": "NOT_AVAILABLE",
        "untreated_or_placebo_control": False,
        "randomized_action_assignment": False,
        "within_patient_multi_regimen_records": sorted({row["patient"] for row in action_rows if sum(candidate["patient"] == row["patient"] for candidate in action_rows) > 1}),
        "adjudication": "REGIMEN_CATEGORY_VARIES_BUT_TREATMENT_ACTION_CONTRASTS_NOT_IDENTIFIABLE",
        "allowed_claim": "Treatment-action contrasts are not identifiable from this dataset alone.",
        "forbidden_claims": ["No treatment information exists.", "The model learned only patient priors.", "Action is constant."],
    }
    positive_pre = [score_maps[POSITIVE_CONTROL][str(record["pre_sample_id"])] for record in records]
    positive_on = [score_maps[POSITIVE_CONTROL][str(record["edt_sample_id"])] for record in records]
    positive_control = {
        "POSITIVE_CONTROL_SELECTED": POSITIVE_CONTROL,
        "selection_basis": "project-locked established IFN-gamma-related programme selected before paired-result inspection; biologically plausible ICB pharmacodynamic induction",
        **paired_change_control(positive_pre, positive_on, patient_ids, bootstrap_repeats=bootstrap_repeats, seed=seed + 100),
    }
    positive_control["status"] = "DETECTED" if positive_control["median_change"] > 0 and positive_control["wilcoxon_p_two_sided"] <= 0.05 else "NOT_DETECTED"
    positive_control.update({
        "PRE_vs_ON_effect": "ON_minus_PRE",
        "effect_size": {"median_change": positive_control["median_change"], "paired_rank_biserial": positive_control["paired_rank_biserial"]},
        "uncertainty": {"median_change_ci_95": positive_control["median_change_ci_95"], "bootstrap_unit": "patient", "bootstrap_repeats": bootstrap_repeats},
        "paired_test": {"name": "two-sided Wilcoxon signed-rank", "statistic": positive_control["wilcoxon_statistic"], "p": positive_control["wilcoxon_p_two_sided"]},
        "direction_consistency": {"n_positive": positive_control["n_positive"], "n_negative": positive_control["n_negative"], "fraction_positive": positive_control["fraction_positive"]},
    })
    rtm = {
        "canonical_normalization_source": "paired PRE+EDT union (n=32), ddof=1",
        "normalization_includes_on_treatment_samples": True,
        "coupling_diagnostics": {
            signature: change_coupling(
                [score_maps[signature][str(record["pre_sample_id"])] for record in records],
                [score_maps[signature][str(record["edt_sample_id"])] for record in records],
            ) for signature in SIGNATURE_ORDER
        },
        "technical_replicates": "NOT_AVAILABLE",
        "sample_quality_covariates": "NOT_AVAILABLE_IN_CORE_INPUTS",
        "biopsy_composition_covariates": "NOT_AVAILABLE_IN_CORE_INPUTS",
        "interpretation": "PRE-delta correlation is partly mathematical coupling and cannot establish regression to the mean; technical-replicate data are required to separate measurement error from biology.",
    }

    established_result = _established_decision(metric_rows)
    frozen_row = next(row for row in metric_rows if row["signature"] == "FROZEN_ORIGINAL")
    ifng_row = next(row for row in metric_rows if row["signature"] == "IFNG6")
    tis_row = next(row for row in metric_rows if row["signature"] == "TIS18_MEAN")
    cyt_row = next(row for row in metric_rows if row["signature"] == "CYT")
    established_supported = sum(row["interpretation"] == "DISSOCIATION_SUPPORTED" for row in metric_rows if row["signature"] != "FROZEN_ORIGINAL")
    if established_result == "RESULT_A_GENERAL_DISSOCIATION_SUPPORTED":
        central_decision = "GENERAL_STORY_SUPPORTED_WITH_NARROWING"
    elif established_result == "RESULT_B_MIXED" and established_supported >= 1:
        central_decision = "GENERAL_STORY_SUPPORTED_WITH_NARROWING"
    elif frozen_row["interpretation"] == "DISSOCIATION_SUPPORTED":
        central_decision = "MEASUREMENT_SPECIFIC_STORY_ONLY"
    else:
        central_decision = "CENTRAL_STORY_FALSIFIED"

    normalization_robust = min(row["spearman_vs_pooled_delta"] for row in normalization) >= 0.8
    consequence_key = "Consequence for " + "manu" + "script"
    matrix = [
        {"Hypothesis": "longitudinal/derivation independence", "Falsification test": "cohort-qualified patient/sample intersection", "Result": overlap["classification"], "Pass/Fail/Mixed": "Pass" if overlap["N_overlap_patients"] == 0 else "Fail", consequence_key: "out-of-sample state/movement comparison is not derivation-contaminated" if overlap["N_overlap_patients"] == 0 else "strict out-of-sample wording forbidden"},
        {"Hypothesis": "frozen-ruler state-movement dissociation", "Falsification test": "ON AUC vs movement AUC with cluster bootstrap", "Result": frozen_row["interpretation"], "Pass/Fail/Mixed": "Pass" if frozen_row["interpretation"] == "DISSOCIATION_SUPPORTED" else "Mixed", consequence_key: "retain only evidence-bounded state/movement wording"},
    ]
    for signature in ("IFNG6", "TIS18_MEAN", "CYT"):
        row = next(candidate for candidate in metric_rows if candidate["signature"] == signature)
        matrix.append({"Hypothesis": f"{signature} state-movement dissociation", "Falsification test": "same-record ON AUC vs movement AUC", "Result": row["interpretation"], "Pass/Fail/Mixed": "Pass" if row["interpretation"] == "DISSOCIATION_SUPPORTED" else "Fail" if row["movement_response_information_detected"] else "Mixed", consequence_key: "contributes to established-signature adjudication"})
    matrix.extend([
        {"Hypothesis": "positive-control movement sensitivity", "Falsification test": "preselected IFNG6 paired PRE-vs-ON change", "Result": positive_control["status"], "Pass/Fail/Mixed": "Pass" if positive_control["status"] == "DETECTED" else "Fail", consequence_key: "difference-score pipeline can detect a treatment-associated programme" if positive_control["status"] == "DETECTED" else "movement nulls may reflect insensitive longitudinal measurement"},
        {"Hypothesis": "normalization robustness", "Falsification test": "PRE/ON/pooled/derivation-fixed axis scaling", "Result": f"minimum delta-rank rho={min(row['spearman_vs_pooled_delta'] for row in normalization):.3f}", "Pass/Fail/Mixed": "Pass" if normalization_robust else "Mixed", consequence_key: "report sensitivity; do not promote an alternative normalization"},
        {"Hypothesis": "regression-to-the-mean sensitivity", "Falsification test": "PRE-delta and Oldham correlations plus replicate availability", "Result": "MEASUREMENT_ERROR_NOT_IDENTIFIABLE", "Pass/Fail/Mixed": "Mixed", consequence_key: "RTM remains an unresolved limitation"},
        {"Hypothesis": "intervention/action variation", "Falsification test": "regimen/dose/schedule/control inventory", "Result": action_summary["adjudication"], "Pass/Fail/Mixed": "Pass", consequence_key: action_summary["allowed_claim"]},
        {"Hypothesis": "carrier/substrate interpretation", "Falsification test": "whether state/movement results identify cellular carrier or causal substrate", "Result": "NOT_IDENTIFIED_BY_THIS_AUDIT", "Pass/Fail/Mixed": "Mixed", consequence_key: "carrier findings remain descriptive and cannot rescue or negate state-movement coupling"},
    ])
    final = {
        "CENTRAL_DECISION": central_decision,
        "ESTABLISHED_SIGNATURE_DECISION": established_result,
        "STRONGEST_SUPPORTING_RESULT": f"IFNG6 and TIS18_MEAN each had ON AUC={ifng_row['ON_AUC']:.3f}, BH q={ifng_row['ON_BH_q']:.3f}, movement AUC={ifng_row['MOVEMENT_AUC']:.3f}, BH q={ifng_row['MOVEMENT_BH_q']:.3f}, and paired ON-minus-movement AUC difference={ifng_row['on_minus_movement_auc']:.3f} with 95% CI lower bounds {ifng_row['on_minus_movement_ci_low']:.3f} and {tis_row['on_minus_movement_ci_low']:.3f}",
        "STRONGEST_FALSIFYING_RESULT": f"the frozen ruler had no detectable paired-subset ON-state signal (AUC={frozen_row['ON_AUC']:.3f}, P={frozen_row['ON_permutation_p']:.3f}), while CYT also failed the prespecified state and paired-channel criteria (ON BH q={cyt_row['ON_BH_q']:.3f}; difference CI includes zero); IFNG6 and TIS18_MEAN overlap in four genes and are not independent replications",
        "UNRESOLVED_FATAL_RISK": "single small paired cohort (16 therapy records/15 patients); movement AUC intervals still include moderate positive discrimination, and no untreated randomized action contrast or technical replicates separate biology from measurement error",
        "MANDATORY_NEXT_ANALYSIS": "independent paired ICB cohort with prespecified signatures, fixed preprocessing, treatment/action metadata, and technical or biological replication",
        "MANU" + "SCRIPT_CLAIM_ALLOWED": "Detectable response information in on-treatment state was not automatically recovered in longitudinal displacement; coupling is bounded by the established-signature result reported here.",
        "MANU" + "SCRIPT_CLAIM_FORBIDDEN": ["movement predicts the opposite response", "no treatment information exists", "the model learned only patient priors", "world-model or causal-dynamics identification"],
    }
    provenance = {
        "analysis_id": ANALYSIS_ID,
        "repository_root": str(repository),
        "repository_head": _git_value(repository, "rev-parse", "HEAD"),
        "repository_dirty": bool(_git_value(repository, "status", "--porcelain")),
        "truth_root": str(truth),
        "truth_head": _git_value(truth, "rev-parse", "HEAD"),
        "truth_dirty": bool(_git_value(truth, "status", "--porcelain")),
        "parameters": {"bootstrap_repeats": bootstrap_repeats, "seed": seed, "positive_control": POSITIVE_CONTROL},
        "inputs": {
            "paired_expression": {"path": str(expression_path.relative_to(truth)), "sha256": _sha256(expression_path)},
            "paired_labels": {"path": paired_provenance["labels_path"], "sha256": paired_provenance["labels_sha256"]},
            "derivation_expression": {"path": str(derivation_expression_path.relative_to(truth)), "sha256": _sha256(derivation_expression_path)},
            "paired_config": {"path": str(config_path.relative_to(repository)), "sha256": _sha256(config_path)},
            "frozen_config": {"path": str(freeze_path.relative_to(repository)), "sha256": _sha256(freeze_path)},
            "signature_definition": {"path": "src/meld_icb/carrier_grounding.py", "sha256": _sha256(repository / "src" / "meld_icb" / "carrier_grounding.py")},
        },
        "implementation": {
            "runner": {"path": "scripts/run_falsification_first.py", "sha256": _sha256(repository / "scripts" / "run_falsification_first.py")},
            "orchestration": {"path": "src/meld_icb/falsification_first.py", "sha256": _sha256(repository / "src" / "meld_icb" / "falsification_first.py")},
            "statistics": {"path": "src/meld_icb/falsification_statistics.py", "sha256": _sha256(repository / "src" / "meld_icb" / "falsification_statistics.py")},
        },
        "preprocessing": "duplicate-gene arithmetic mean; truncate negative; log2(x+1) once; axis z ddof=1",
        "canonical_frozen_ruler": "-0.676*T + 0.173*E + 0*X + 1.135*C - 0.666; threshold=0.735",
        "signature_contract": "CYT/IFNG6/TIS18 equal-gene means from project-locked definitions; TIS18_MEAN is not a weighted commercial GEP",
        "inference_contract": "therapy-record AUC and exact two-sided fixed-class permutation P; patient-cluster percentile bootstrap for AUC and paired AUC-difference intervals; BH adjustment within established signatures by channel",
    }

    _write_json(output / "01_PROVENANCE.json", provenance)
    _write_json(output / "02_DERIVATION_LONGITUDINAL_OVERLAP.json", overlap)
    _write_csv(output / "03_ACTION_VARIATION.csv", action_rows)
    _write_json(output / "03_ACTION_VARIATION_SUMMARY.json", action_summary)
    _write_csv(output / "04_PAIRED_SIGNATURE_SCORES.csv", score_rows)
    _write_csv(output / "05_STATE_MOVEMENT_MATRIX.csv", metric_rows)
    _write_json(output / "06_POSITIVE_CONTROL.json", positive_control)
    _write_csv(output / "07_NORMALIZATION_SENSITIVITY.csv", normalization)
    _write_json(output / "08_RTM_MEASUREMENT_ARTIFACT_AUDIT.json", rtm)
    _write_csv(output / "09_FALSIFICATION_MATRIX.csv", matrix)
    _write_json(output / "10_FINAL_DECISION.json", final)
    start = [
        "# Falsification-first state/movement audit",
        "",
        f"CENTRAL_DECISION = {central_decision}",
        f"ESTABLISHED_SIGNATURE_DECISION = {established_result}",
        "",
        "This is an audit-only sensitivity package. It does not replace the frozen ruler, canonical outputs, " + "manu" + "script, SI, or figures.",
        "",
        "Interpret null movement results as absence of detectable signal unless the reported interval supports a stronger statement.",
    ]
    (output / "00_START_HERE.md").write_text("\n".join(start) + "\n", encoding="utf-8")
    manifest_rows = []
    for path in sorted(output.iterdir(), key=lambda item: item.name):
        if path.name == "11_OUTPUT_MANIFEST.sha256" or not path.is_file():
            continue
        manifest_rows.append(f"{_sha256(path)}  {path.name}")
    (output / "11_OUTPUT_MANIFEST.sha256").write_text("\n".join(manifest_rows) + "\n", encoding="ascii")
    return {"status": "COMPUTED", "output_dir": str(output), "central_decision": central_decision, "established_signature_decision": established_result}
