from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from meld_icb.object_adapters import all_object_ids, build_adapter_payload, load_object_map, run_object_adapter
from scripts.public_runner import produce_computed


ROOT = Path(__file__).resolve().parents[1]


def _write(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload), encoding="utf-8")


def _valid_b() -> dict:
    return {
        "status": "COMPUTED",
        "b2_movement": {"n_records": 16, "n_patients": 15, "median_record_level_cosine": 0.6993631612147504},
        "b3_response_specificity": {"delta_score_auc": 0.3272727272727272},
        "b4_incremental_foldsafe": {"P_two_sided": 0.03094482421875},
    }


def _valid_c() -> dict:
    return {"verdict": "GATES_PASSED", "n_cells": 10363, "gates": {"A_lymphoid_share": {"observed": 0.8629283905029297}, "B_malignant_share": {"observed": 0.0007919137715362012}, "C_top_two_carriers": {"observed": ["CD8_T", "CD4_T"]}}}


def test_map_and_manifest_cover_exactly_25_objects():
    assert len(all_object_ids(ROOT)) == 25
    mapping = load_object_map(ROOT)
    with (ROOT / "manifests" / "producer_status.tsv").open(encoding="utf-8") as handle:
        statuses = {row["object_id"]: row["status"] for row in csv.DictReader(handle, delimiter="\t")}
    assert set(statuses) == set(mapping)
    assert all(statuses[key] == mapping[key]["default_status"] for key in mapping)
    with (ROOT / "manifests" / "producer_status.tsv").open(encoding="utf-8") as handle:
        rows = {row["object_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    assert rows["ST03"]["claim_ids"] == "C2;C3;C4"
    assert rows["ST03"]["required_cohorts"] == "PRJEB23709"
    assert rows["ST04"]["claim_ids"] == "C5"
    assert rows["ST04"]["required_cohorts"] == "Sade-Feldman"


@pytest.mark.parametrize("object_id", all_object_ids(ROOT))
def test_every_object_has_adapter_status(object_id, tmp_path):
    payload = build_adapter_payload(object_id, shared_dir=tmp_path, repo_root=ROOT)
    assert payload["object_id"] == object_id
    assert payload["status"]


def test_missing_mapped_chain_is_hold_and_exit_two(tmp_path):
    assert run_object_adapter("Figure3", shared_dir=tmp_path, output_dir=tmp_path, repo_root=ROOT) == 2
    payload = json.loads((tmp_path / "Figure3.status.json").read_text(encoding="utf-8"))
    assert payload["status"] == "HOLD_COMPONENT_NOT_AVAILABLE"
    assert payload["chain_status"] == "HOLD"


def test_expected_value_mismatch_is_hold(tmp_path):
    _write(tmp_path / "B_PUBLIC_NUMERIC_OUTPUT.json", {**_valid_b(), "b2_movement": {"n_records": 99, "n_patients": 15, "median_record_level_cosine": 0.6993631612147504}})
    payload = build_adapter_payload("Figure3", shared_dir=tmp_path, repo_root=ROOT)
    assert payload["status"] == "HOLD_EXPECTED_VALUE_MISMATCH"
    assert payload["expected_validation"]["details"]["b2_n_records"]["tolerance"] == 1e-12


def test_valid_b_and_c_aggregate_projections(tmp_path):
    _write(tmp_path / "B_PUBLIC_NUMERIC_OUTPUT.json", _valid_b())
    _write(tmp_path / "C_PUBLIC_NUMERIC_OUTPUT.json", _valid_c())
    assert build_adapter_payload("Figure3", shared_dir=tmp_path, repo_root=ROOT)["status"] == "COMPUTED"
    assert build_adapter_payload("Figure4", shared_dir=tmp_path, repo_root=ROOT)["status"] == "COMPUTED"


def test_a_aggregate_projection_validates_cohorts_pooled_and_identity(tmp_path):
    _write(tmp_path / "A_PUBLIC_NUMERIC_OUTPUT.json", {"status": "COMPUTED", "cohorts": {"PRJEB23709": {"status": "COMPUTED", "metrics": {"auc": 0.8441558441558441}}, "MORRISON-1-public": {"status": "COMPUTED", "metrics": {"auc": 0.825}}, "MGH": {"status": "COMPUTED", "metrics": {"auc": 0.675}}}, "pooled": {"n_records": 85, "auc": 0.7705345501955672, "ci_95": [0.6538461538461539, 0.8800521512385919], "p_two_sided": 0.000999000999000999}, "identity": {"n_records": 85, "sha256": "bedb09abc0ee1338fc904f834cda957751a79202cff9102f1c5fe986a11ddd11"}})
    payload = build_adapter_payload("Figure2", shared_dir=tmp_path, repo_root=ROOT)
    assert payload["status"] == "COMPUTED"
    assert payload["expected_validation"]["status"] == "PASS"
    assert payload["expected_validation"]["details"]["pooled_p_two_sided"]["tolerance"] == 1e-12
    assert payload["expected_validation"]["details"]["identity_sha256"]["error"] == 0


def test_d_partial_projection_validates_gse78220_without_upgrading_verdict(tmp_path):
    _write(tmp_path / "D_PUBLIC_NUMERIC_OUTPUT.json", {"status": "HOLD", "cohorts": {"GSE78220": {"status": "COMPUTED", "populations": {"primary": {"metrics": {"auc": 0.47619047619047616, "ci_95": [0.24404761904761904, 0.7083333333333334], "p_two_sided": 0.8573142685731426, "n_patients": 26}}}}}})
    payload = build_adapter_payload("Figure5", shared_dir=tmp_path, repo_root=ROOT)
    assert payload["status"] == "COMPUTED_PRIMARY_PLUS_COMPARISON_ONLY"
    assert payload["expected_validation"]["status"] == "PASS"


def test_d_partial_projection_mismatch_holds_expected_value(tmp_path):
    _write(tmp_path / "D_PUBLIC_NUMERIC_OUTPUT.json", {"status": "HOLD", "cohorts": {"GSE78220": {"status": "COMPUTED", "populations": {"primary": {"metrics": {"auc": 0.1, "ci_95": [0.2, 0.3], "p_two_sided": 0.4, "n_patients": 26}}}}}})
    payload = build_adapter_payload("Figure5", shared_dir=tmp_path, repo_root=ROOT)
    assert payload["status"] == "HOLD_EXPECTED_VALUE_MISMATCH"


def test_d_mixed_projection_mismatch_holds_even_when_top_level_computed(tmp_path):
    _write(tmp_path / "D_PUBLIC_NUMERIC_OUTPUT.json", {"status": "COMPUTED", "cohorts": {"GSE78220": {"status": "COMPUTED", "populations": {"primary": {"metrics": {"auc": 0.1, "ci_95": [0.2, 0.3], "p_two_sided": 0.4, "n_patients": 26}}}}}})
    payload = build_adapter_payload("Figure5", shared_dir=tmp_path, repo_root=ROOT)
    assert payload["status"] == "HOLD_EXPECTED_VALUE_MISMATCH"


def test_st07_is_retired_and_st08_policy_contract(tmp_path):
    payload = build_adapter_payload("ST07", shared_dir=tmp_path, repo_root=ROOT)
    assert payload["status"] == "RETIRED"
    assert payload["active_projection"] is False
    assert "superseded" in payload["retirement_reason"]
    excluded = build_adapter_payload("ST08", shared_dir=tmp_path, repo_root=ROOT)
    assert excluded["status"] == "EXCLUDED_BY_POLICY"
    assert excluded["schema_only"]["omitted_patient_details_included"] is False
    assert excluded["schema_only"]["internal_ids_included"] is False
    assert "patient_id" not in json.dumps(excluded["schema_only"])


def test_figure1_conceptual_contract_has_no_hold_reason(tmp_path):
    payload = build_adapter_payload("Figure1", shared_dir=tmp_path, repo_root=ROOT)
    assert payload["status"] == "CONCEPTUAL_ASSET_CONTRACT"
    assert "hold_reason" not in payload


def test_old_raw_bypass_is_fail_closed(tmp_path):
    assert produce_computed("Figure3", [], [], "figure3", ["raw.tsv"], output_dir=str(tmp_path)) == 2
    payload = json.loads((tmp_path / "Figure3.status.json").read_text(encoding="utf-8"))
    assert payload["status"] == "HOLD_RAW_RECOMPUTE_RETIRED"


def test_wrappers_do_not_import_core_producers_or_statistics():
    scripts = [path for path in (ROOT / "scripts").glob("reproduce_*.py") if path.name != "reproduce_all.py"]
    assert len(scripts) == 25
    for script in scripts:
        text = script.read_text(encoding="utf-8")
        assert "meld_icb.producers" not in text
        assert "statistics" not in text
