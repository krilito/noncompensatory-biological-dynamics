"""Scientific-contract tests. These are not existence or smoke checks."""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from meld_icb.cli import main as cli_main
from meld_icb.config import load_freeze
from meld_icb.incremental_value import fold_safe_incremental_metrics, influence_leave_one_patient_reruns
from meld_icb.paired_movement import paired_movement_metrics
from meld_icb.response_specificity import response_specificity_metrics


ROOT = Path(__file__).resolve().parents[1]


def _csv(name: str) -> list[dict[str, str]]:
    with (ROOT / "source_data" / name).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _table6() -> dict[str, dict[str, str]]:
    rows = {}
    with (ROOT / "source_data" / "supplementary_tables" / "table6.csv").open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            rows[row["Audit"].split()[0]] = row
    return rows


def test_claim_vector_is_noncompensatory_and_negative_findings_remain():
    table = _table6()
    assert "CONFLICT" in table["B0"]["Decision"]
    assert table["B1"]["Decision"] == "SUPPORTED"
    assert table["B2"]["Decision"] == "SUPPORTED"
    assert table["B3"]["Decision"] == "NOT SUPPORTED"
    assert table["B4"]["Decision"] == "NOT SUPPORTED"
    assert table["B5"]["Decision"] == "DESCRIPTIVE"
    assert table["Treatment"]["Decision"] == "NOT TESTED"
    mgh = next(row for row in _csv("figure2_auc.csv") if row["cohort"] == "MGH sample-level")
    assert mgh["role"] == "supportive"
    assert float(mgh["auc_canonical"]) < 0.7


def test_b0_historical_lineage_is_not_washed_out():
    rows = {row["cohort"]: row for row in _csv("figure2_auc.csv")}
    assert float(rows["PRJEB23709"]["auc_canonical"]) != float(rows["PRJEB23709"]["auc_historical"])
    assert float(rows["MGH sample-level"]["auc_canonical"]) != float(rows["MGH sample-level"]["auc_historical"])
    assert float(rows["MGH patient-level"]["p_canonical"]) != float(rows["MGH patient-level"]["p_historical"])


def test_b2_record_and_cohort_cosines_are_distinct_objects():
    summary = _csv("figure3_b2_summary.csv")[0]
    record = float(summary["median_record_level_cosine"])
    cohort = float(summary["cohort_transition_cosine"])
    assert record != cohort
    assert 0.6 < record < 0.8
    assert 0.98 < cohort <= 1.0
    freeze = load_freeze(ROOT)
    records = [
        {"patient_id": 1, "pre_sample_id": "pre1", "edt_sample_id": "edt1"},
        {"patient_id": 2, "pre_sample_id": "pre2", "edt_sample_id": "edt2"},
    ]
    z_axes = {
        "pre1": {"T": 0.0, "E": 0.0, "X": 0.0, "C": 0.0},
        "edt1": {"T": 0.0, "E": 1.0, "X": 0.0, "C": 0.0},
        "pre2": {"T": 0.0, "E": 0.0, "X": 0.0, "C": 0.0},
        "edt2": {"T": 1.0, "E": 0.0, "X": 0.0, "C": 0.0},
    }
    result = paired_movement_metrics(records, z_axes, {"T": 0.0, "E": 1.0, "X": 0.0, "C": 0.0}, freeze)
    assert "median_record_level_cosine" in result
    assert "cohort_transition_cosine_to_GSE91061" in result
    assert result["median_record_level_cosine"] != result["cohort_transition_cosine_to_GSE91061"]


def test_b3_does_not_invert_unfavorable_auc():
    records = [{"patient_id": index, "y_true": int(index >= 3)} for index in range(6)]
    rows = [{"delta_boundary": float(5 - index), "cosine_to_reference": float(5 - index)} for index in range(6)]
    result = response_specificity_metrics(records, rows, bootstrap_repeats=20, permutation_repeats=20)
    assert result["delta_score_auc"] < 0.5
    assert result["canonical_status"] == "NOT_SUPPORTED"
    committed = next(row for row in _csv("figure3_b3.csv") if row["readout"] == "boundary change")
    assert float(committed["auc"]) < 0.5
    assert int(committed["n_R"]) == 11
    assert int(committed["n_NR"]) == 5


def test_b4_fold_safe_excludes_held_out_patient_from_normalization():
    records = []
    raw_axes: dict[str, dict[str, float]] = {}
    for index in range(16):
        patient = index if index < 15 else 0
        pre = f"pre{index}"
        edt = f"edt{index}"
        records.append(
            {
                "patient_id": patient,
                "pre_sample_id": pre,
                "edt_sample_id": edt,
                "y_true": int(index % 2),
            }
        )
        raw_axes[pre] = {"T": float(index), "E": float(index + 1), "X": float(index + 2), "C": float(index + 3)}
        raw_axes[edt] = {"T": float(index + 4), "E": float(index + 5), "X": float(index + 6), "C": float(index + 7)}
    result = fold_safe_incremental_metrics(records, raw_axes, load_freeze(ROOT), canonical_scope=True)
    assert result["axis_normalization_train_fold_only"] is True
    assert result["n_patients"] == 15
    assert result["canonical_status"] in {"SUPPORTED", "NOT_SUPPORTED"}
    assert result["loss_sign_convention"].startswith("positive delta_L favors augmented")


def test_b4_influence_rebuilds_instead_of_deleting_finished_losses():
    records = []
    raw_axes: dict[str, dict[str, float]] = {}
    for index in range(8):
        pre = f"pre{index}"
        edt = f"edt{index}"
        records.append({"patient_id": index, "pre_sample_id": pre, "edt_sample_id": edt, "y_true": int(index < 4)})
        raw_axes[pre] = {"T": 0.2 * index, "E": 1.0 - 0.1 * index, "X": 0.05 * index, "C": 0.3 * index}
        raw_axes[edt] = {"T": 0.1 * index, "E": 1.2 - 0.05 * index, "X": 0.08 * index, "C": 0.35 * index}
    complete = fold_safe_incremental_metrics(records, raw_axes, load_freeze(ROOT), canonical_scope=False)
    influence = influence_leave_one_patient_reruns(records, raw_axes, load_freeze(ROOT))
    assert influence["rebuilds_folds"] is True
    assert influence["n_reruns"] == 8
    omitted = influence["rows"][0]["omitted_patient"]
    naive = sum(item["delta_L"] for item in complete["patient_delta_L"] if item["patient_id"] != omitted) / 7
    rebuilt = next(row["mean_delta_L"] for row in influence["rows"] if row["omitted_patient"] == omitted)
    assert abs(naive - rebuilt) > 1e-12

    patients = _csv("figure3_b4_patients.csv")
    committed = {int(row["omitted_patient"]): float(row["mean_delta_L"]) for row in _csv("figure3_influence.csv")}
    for omitted_id in (17, 13):
        remaining = [float(row["patient_delta_L_X0_minus_X1"]) for row in patients if int(row["patient_id"]) != omitted_id]
        naive_committed = sum(remaining) / len(remaining)
        assert abs(naive_committed - committed[omitted_id]) > 1e-4


def test_b4_committed_sign_is_negative_and_not_supported():
    summary = _csv("figure3_b4_summary.csv")[0]
    assert float(summary["mean_delta_L"]) < 0
    assert float(summary["L_augmented"]) > float(summary["L_baseline"])
    assert int(summary["n_sign_patterns"]) == 32768
    assert {int(row["patient_id"]) for row in _csv("figure3_b4_patients.csv")} == {
        1, 8, 10, 13, 17, 19, 25, 27, 29, 33, 35, 36, 41, 47, 48
    }
    assert sum(int(row["n_heldout_records"]) for row in _csv("figure3_b4_patients.csv")) == 16


def test_b5_t_carrier_shares_remain_descriptive():
    rows = [row for row in _csv("figure4.csv") if row["axis"] == "T"]
    shares = {row["lineage"]: float(row["axis_signal_share"]) for row in rows}
    lymphoid = sum(shares[name] for name in ("CD4_T", "CD8_T", "NK", "B_cell", "Treg", "Plasma"))
    assert abs(lymphoid - 0.86293) < 5e-6
    assert abs(shares["Malignant"] - 0.0007919137715362) < 1e-12


def test_b6_roles_are_not_upgraded_to_external_validation():
    rows = {row["dataset"]: row for row in _csv("figure5.csv")}
    assert rows["GSE78220"]["evidence_role"] == "NEAR_NULL"
    assert "COMPARISON" in rows["IMvigor210"]["evidence_role"] or "UNRESOLVED" in rows["IMvigor210"]["evidence_role"]
    assert "CARRIER_CONTEXT" in rows["Sade-Feldman"]["evidence_role"]


def test_imvigor_sample_table_is_not_redistributed():
    assert not (ROOT / "source_data" / "figure5_imvigor_samples.csv").exists()
    plotter = (ROOT / "scripts" / "plot_figure5.py").read_text(encoding="utf-8")
    assert "figure5_imvigor_samples" not in plotter
    assert "ACCESSION_ONLY" in (ROOT / "manifests" / "cohorts.tsv").read_text(encoding="utf-8")


def test_citation_identifiers_are_current_and_not_placeholders():
    # After Zenodo publishes v1.0.1, cite that exact release and retain no placeholders.
    text = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
    assert "arxiv:TBD" not in text
    assert "identifiers:" not in text
    assert 'version: "1.0.1"' in text
    assert 'doi: "10.5281/zenodo.22273685"' in text
    assert "date-released: 2026-09-03" in text


def test_cli_loads_freeze_without_private_representation_key(capsys):
    assert cli_main(["--repo-root", str(ROOT)]) == 0
    captured = capsys.readouterr().out
    assert "A_FROZEN_STATE_TRANSFER_V1" in captured
    assert "0.735" in captured
