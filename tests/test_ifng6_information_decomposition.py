from __future__ import annotations

import math

import pytest

from meld_icb.falsification_statistics import clustered_auc_analysis, clustered_correlation_analysis
from meld_icb.ifng6_decomposition_models import lopo_incremental_information, residualized_change
from meld_icb.ifng6_information_decomposition import (
    _assert_locked_regression,
    _decomposition_decision,
    _scope_contract,
)


def test_clustered_auc_analysis_uses_one_draw_for_paired_differences() -> None:
    labels = [0, 0, 0, 1, 1, 1]
    patients = ["a", "b", "c", "d", "e", "e"]
    channels = {
        "PRE": [0.0, 0.1, 0.2, 0.7, 0.8, 0.9],
        "ON": [0.0, 0.2, 0.1, 0.8, 0.9, 1.0],
    }
    summary, draws = clustered_auc_analysis(
        channels, labels, patients, [("ON", "PRE")], bootstrap_repeats=200, seed=4
    )
    assert summary["n_patients"] == 5
    assert draws
    assert all(math.isclose(row["ON_MINUS_PRE"], row["AUC_ON"] - row["AUC_PRE"]) for row in draws)


def test_decomposition_can_forbid_record_label_permutation() -> None:
    summary, _ = clustered_auc_analysis(
        {"PRE": [0.0, 0.1, 0.8, 0.9]},
        [0, 0, 1, 1],
        ["a", "b", "c", "d"],
        [],
        bootstrap_repeats=50,
        seed=5,
        compute_record_label_permutation=False,
    )
    assert summary["channels"]["PRE"]["auc_permutation_p_two_sided"] is None
    assert summary["permutation_unit_note"].startswith("not run")


def test_residualization_is_label_free_and_lopo_holds_out_patient_cluster() -> None:
    pre = [0.0, 1.0, 2.0, 3.0, 3.5]
    on = [0.2, 1.4, 2.1, 3.7, 4.0]
    patients = ["1", "2", "3", "4", "4"]
    result = residualized_change(pre, on, patients)
    assert len(result["full_sample"]["residuals"]) == 5
    held_four = next(fold for fold in result["lopo"]["folds"] if fold["held_out_patient"] == "4")
    assert held_four["n_held_out_records"] == 2
    assert abs(clustered_correlation_analysis(
        pre, result["full_sample"]["residuals"], patients, bootstrap_repeats=50, seed=2
    )["pearson_r"]) < 1e-12


def _decision_fixture(*, pre: tuple[float, float], on: tuple[float, float], raw: tuple[float, float], residual: tuple[float, float]):
    channels = {
        name: {"auc": values[0], "auc_ci_95": [values[1], 1.0]}
        for name, values in {"PRE": pre, "ON": on, "RAW_DELTA": raw, "RESIDUAL_LOPO": residual}.items()
    }
    aucs = {name: values["auc"] for name, values in channels.items()}
    comparisons = {
        f"{left}_MINUS_{right}": {"estimate": aucs[left] - aucs[right], "ci_95": [-0.1, 0.5]}
        for left in ("PRE", "ON") for right in ("RAW_DELTA", "RESIDUAL_LOPO")
    }
    return {"channels": channels, "comparisons": comparisons}


def test_interpretation_rules_are_frozen_before_results() -> None:
    assert "Decision precedence" in _scope_contract()
    outcome_a = _decision_fixture(pre=(0.84, 0.57), on=(0.91, 0.70), raw=(0.42, 0.15), residual=(0.60, 0.30))
    assert _decomposition_decision(outcome_a) == "STATE_INFORMATION_PERSISTS_CHANGE_NOT_READOUT"
    outcome_c = _decision_fixture(pre=(0.84, 0.57), on=(0.91, 0.70), raw=(0.42, 0.15), residual=(0.80, 0.56))
    assert _decomposition_decision(outcome_c) == "RAW_CHANGE_ARTIFACT_OR_PARAMETERIZATION_LIMITATION"


def test_locked_ifng6_regression_gate_rejects_numeric_drift() -> None:
    primary = {
        "channels": {
            name: {"auc": auc, "auc_ci_95": [low, high], "auc_permutation_p_two_sided": p}
            for name, auc, low, high, p in [
                ("PRE", 0.8, 0.6, 1.0, 0.04),
                ("ON", 0.9, 0.7, 1.0, 0.01),
                ("RAW_DELTA", 0.4, 0.1, 0.7, 0.6),
            ]
        },
        "comparisons": {"ON_MINUS_RAW_DELTA": {"estimate": 0.5, "ci_95": [0.1, 0.8]}},
    }
    locked = {
        "PRE_AUC": "0.8", "PRE_CI_low": "0.6", "PRE_CI_high": "1.0", "PRE_permutation_p": "0.04",
        "ON_AUC": "0.9", "ON_CI_low": "0.7", "ON_CI_high": "1.0", "ON_permutation_p": "0.01",
        "MOVEMENT_AUC": "0.4", "MOVEMENT_CI_low": "0.1", "MOVEMENT_CI_high": "0.7", "MOVEMENT_permutation_p": "0.6",
        "on_minus_movement_auc": "0.5", "on_minus_movement_ci_low": "0.1", "on_minus_movement_ci_high": "0.8",
    }
    _assert_locked_regression(primary, locked)
    locked["ON_AUC"] = "0.89"
    with pytest.raises(ValueError, match="ON_AUC"):
        _assert_locked_regression(primary, locked)


def test_incremental_models_are_patient_lopo_and_return_uncertainty() -> None:
    patient_ids = [str(value) for value in range(1, 16)] + ["13"]
    labels = [0, 0, 0, 0, 0] + [1] * 11
    rows = []
    for index, (patient, label) in enumerate(zip(patient_ids, labels)):
        pre = float(index) / 4.0
        on = pre + (0.2 if label else 0.8)
        rows.append({
            "patient_id": patient,
            "therapy_record_id": f"r{index}",
            "response_binary": label,
            "IFNG6_PRE": pre,
            "IFNG6_ON": on,
            "IFNG6_RAW_DELTA": on - pre,
        })
    result = lopo_incremental_information(rows, bootstrap_repeats=100, seed=8)
    assert set(result["models"]) == {"M_PRE", "M_ON", "M_PRE_ON", "M_ON_DELTA"}
    assert result["bootstrap_repeats_valid"] > 0
    assert result["incremental_on_given_pre"]["status"] in {
        "ADDS_HELD_OUT_INFORMATION", "DEGRADES_HELD_OUT_INFORMATION", "INCREMENTAL_ANALYSIS_UNSTABLE"
    }
