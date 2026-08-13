import json
import os
from pathlib import Path

import pytest

from meld_icb.config import load_freeze
from meld_icb.incremental_value import exact_sign_flip_distribution, fold_safe_incremental_metrics
from meld_icb.paired_movement import cosine, paired_movement_metrics
from meld_icb.response_specificity import response_specificity_metrics
from scripts.run_paired_dynamics import run


ROOT = Path(__file__).resolve().parents[1]


def _record(index: int, patient: int | None = None) -> dict[str, object]:
    patient_id = patient if patient is not None else index
    return {
        "record_id": f"r{index}",
        "patient_id": patient_id,
        "therapy_type": "PD1",
        "pre_sample_id": f"pre{index}",
        "edt_sample_id": f"edt{index}",
        "y_true": int(index % 2),
    }


def test_pairing_metrics_preserve_wilcoxon_and_cosine_contract():
    records = [_record(1), _record(2)]
    z_axes = {
        "pre1": {"T": 0.0, "E": 0.0, "X": 0.0, "C": 0.0},
        "edt1": {"T": 0.0, "E": 1.0, "X": 0.0, "C": 0.0},
        "pre2": {"T": 0.0, "E": 0.0, "X": 0.0, "C": 0.0},
        "edt2": {"T": 0.0, "E": 0.0, "X": 0.0, "C": 1.0},
    }
    freeze = load_freeze(ROOT)
    result = paired_movement_metrics(records, z_axes, {"T": 0.0, "E": 1.0, "X": 0.0, "C": 0.0}, freeze)
    assert result["n_records"] == 2
    assert result["fraction_delta_gt0"] == 1.0
    assert result["wilcoxon_pvalue"] == 0.5
    assert cosine({"T": 0, "E": 1, "X": 0, "C": 0}, {"T": 0, "E": 1, "X": 0, "C": 0}) == 1.0


def test_specificity_retains_small_sample_not_supported_semantics():
    records = [_record(1), _record(2), _record(3), _record(4)]
    rows = [{"delta_boundary": float(index), "cosine_to_reference": float(4 - index)} for index in range(4)]
    result = response_specificity_metrics(records, rows, bootstrap_repeats=20, permutation_repeats=20)
    assert result["canonical_status"] == "NOT_SUPPORTED"
    assert "not evidence of biological absence" in result["interpretation"]


def test_exact_sign_flip_is_exhaustive_and_loss_sign_is_explicit():
    result = exact_sign_flip_distribution([1.0, -1.0])
    assert result["n_sign_patterns"] == 4
    assert result["P_improvement"] == 0.75
    assert result["P_worsening"] == 0.75


def test_b4_rejects_duplicate_sample_leakage():
    records = [_record(index, patient=index if index < 15 else 1) for index in range(16)]
    records[-1]["pre_sample_id"] = records[0]["pre_sample_id"]
    raw_axes = {sample_id: {"T": float(index), "E": float(index + 1), "X": float(index + 2), "C": float(index + 3)} for index, sample_id in enumerate({key for record in records for key in (record["pre_sample_id"], record["edt_sample_id"])})}
    with pytest.raises(ValueError, match="sample leakage"):
        fold_safe_incremental_metrics(records, raw_axes, load_freeze(ROOT))


def test_registry_output_cannot_be_used_as_b_producer_input(tmp_path):
    config = json.loads((ROOT / "configs" / "paired_dynamics.yaml").read_text(encoding="utf-8"))
    config["inputs"]["labels"] = "excluded_registry/PAIRED_RECORD_REGISTRY.csv"
    config_path = tmp_path / "paired.json"
    config_path.write_text(json.dumps(config), encoding="utf-8")
    expression = tmp_path / "expr.tsv"
    expression.write_text("Gene\tPD1_1_PRE\tPD1_1_EDT\nMKI67\t1\t2\n", encoding="utf-8")
    registry = tmp_path / "registry.csv"
    registry.write_text("record_id,patient_id,therapy_type,pre_sample_id,edt_sample_id,y_true\nr,1,PD1,PD1_1_PRE,PD1_1_EDT,1\n", encoding="utf-8")
    config["inputs"]["expression"] = "expr.tsv"
    config["inputs"]["labels"] = "registry.csv"
    result = run(truth_root=tmp_path, config_path=config_path, output_dir=tmp_path / "out")
    assert result["status"] == "HOLD"


def test_missing_b_inputs_are_explicit_hold(tmp_path):
    result = run(truth_root=tmp_path / "missing", config_path=ROOT / "configs" / "paired_dynamics.yaml", output_dir=tmp_path / "out")
    assert result["status"] == "HOLD"
