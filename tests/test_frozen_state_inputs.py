import json
from pathlib import Path

from meld_icb.config import load_freeze
from meld_icb.frozen_state_inputs import (
    _collapse_matrix,
    compare_historical_to_canonical,
    discover_gse91061_inputs,
    mgh_earliest_patient_records,
    run_truth_cohorts,
    write_generic_input,
)


ROOT = Path(__file__).resolve().parents[1]


def test_duplicate_collapse_precedes_single_deposited_scale_transform():
    rows = _collapse_matrix(
        {"A": {"s1": [1.0, 3.0]}, "B": {"s1": [0.0]}},
        ["s1"],
        transform="truncate_negative_then_log2p1",
    )
    # mean(1, 3) -> 2, then log2(2 + 1), not mean(log2(1+1), log2(3+1)).
    assert abs(rows[0]["A"] - 1.584962500721156) < 1e-12
    assert rows[0]["B"] == 0.0


def test_mgh_patient_sensitivity_uses_earliest_day_and_stable_source_order():
    records = [
        {"sample_id": "late", "patient_id": "p1", "passon_day": "10", "boundary_score": 2, "response_binary": 1},
        {"sample_id": "early", "patient_id": "p1", "passon_day": "2", "boundary_score": 1, "response_binary": 1},
        {"sample_id": "no_day", "patient_id": "p2", "passon_day": "", "boundary_score": -1, "response_binary": 0},
    ]
    selected = mgh_earliest_patient_records(records)
    assert [(r["patient_id"], r["sample_id"]) for r in selected] == [("p1", "early"), ("p2", "no_day")]


def test_gse91061_discovery_does_not_backfill_comparison_tables(tmp_path):
    # No truth repository is a clean, explicit HOLD; the adapter never treats
    # fig1_lopo_refit/frozen score tables as producer inputs.
    status = discover_gse91061_inputs(tmp_path)
    assert status["status"] == "HOLD"
    assert Path(status["raw_expression"]["path"]).as_posix().endswith("processed/expression_log2_GSE91061_symbol.csv")
    assert status["canonical_ref_fold_inputs"]["exists"] is False


def test_generic_input_is_audit_local_and_contains_only_contract_fields(tmp_path):
    freeze = load_freeze(ROOT)
    path = write_generic_input(
        [{"sample_id": "s1", "patient_id": "p1", "response": 1, "MKI67": 2.0, "unexpected": "drop"}],
        tmp_path / "rows.tsv",
        freeze=freeze,
    )
    text = path.read_text(encoding="utf-8")
    assert "sample_id" in text and "MKI67" in text and "unexpected" not in text


def _comparison_fixture(path: Path, *, extra_id: str | None = None):
    rows = [
        {"_cohort": "PRJEB23709", "sample_id": "s1", "response_binary": 1, "boundary_score": 1.0},
        {"_cohort": "PRJEB23709", "sample_id": "s2", "response_binary": 0, "boundary_score": -1.0},
    ]
    historical = [
        {"cohort": "Gide / PRJEB23709", "sample_id": "s1", "response_binary_responder1": "0", "boundary_score": "0.9"},
        {"cohort": "Gide / PRJEB23709", "sample_id": "s2", "response_binary_responder1": "1", "boundary_score": "-0.8"},
    ]
    if extra_id:
        historical.append({"cohort": "Gide / PRJEB23709", "sample_id": extra_id, "response_binary_responder1": "1", "boundary_score": "0.1"})
    import csv
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(historical[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(historical)
    return rows


def test_historical_comparison_uses_only_legacy_scores_after_canonical_rows_exist(tmp_path):
    records = _comparison_fixture(tmp_path / "historical.csv")
    result = compare_historical_to_canonical(records, tmp_path / "historical.csv")
    assert result["status"] == "COMPUTED"
    assert result["pooled"]["matched_count"] == 2
    # Historical labels are deliberately inverted in the fixture; canonical
    # labels still determine both AUCs, proving the comparison cannot supply
    # producer labels.
    assert result["pooled"]["auc_canonical"] == 1.0
    assert result["pooled"]["orientation_consistent"] is True


def test_historical_comparison_mapping_mismatch_is_fail_closed(tmp_path):
    records = _comparison_fixture(tmp_path / "historical.csv", extra_id="unexpected")
    result = compare_historical_to_canonical(records, tmp_path / "historical.csv")
    assert result["status"] == "HOLD"
    assert result["cohorts"]["Gide / PRJEB23709"]["extra_ids_count"] == 1


def test_missing_historical_object_is_not_supplied_and_does_not_block_producer(tmp_path):
    records = _comparison_fixture(tmp_path / "historical.csv")
    assert compare_historical_to_canonical(records, None)["status"] == "NOT_SUPPLIED"
    result = run_truth_cohorts(
        truth_root=tmp_path / "empty_truth",
        audit_dir=tmp_path / "audit",
        config_path=ROOT / "configs" / "frozen_state_transfer.yaml",
        historical_comparison=tmp_path / "missing.csv",
    )
    assert result["status"] == "HOLD"
    assert result["historical_to_canonical_drift"]["status"] == "HOLD"
