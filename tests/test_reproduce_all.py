import json
import csv
import subprocess
import sys
from pathlib import Path

from scripts.reproduce_all import PRODUCERS
import scripts.reproduce_all as reproduce_all

ROOT = Path(__file__).resolve().parents[1]


def test_analysis_mode_without_data_is_non_success():
    completed = subprocess.run(
        [sys.executable, "scripts/reproduce_all.py", "--mode", "analysis"],
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 2
    assert (ROOT / "figures" / "reproduced" / "Figure1.status.json").exists()


def test_analysis_allow_hold_without_data():
    completed = subprocess.run(
        [sys.executable, "scripts/reproduce_all.py", "--mode", "analysis", "--allow-hold"],
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0
    payload = json.loads((ROOT / "figures" / "reproduced" / "Figure1.status.json").read_text(encoding="utf-8"))
    assert payload["status"] == "CONCEPTUAL_ASSET_CONTRACT"


def test_figures_mode_redraws_without_truth_root():
    completed = subprocess.run(
        [sys.executable, "scripts/reproduce_all.py", "--mode", "figures"],
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0
    assert (ROOT / "figures" / "reproduced" / "Figure2.pdf").exists()
    assert (ROOT / "figures" / "reproduced" / "Figure5.pdf").exists()


def test_fail_on_missing_is_hold_exit_two():
    completed = subprocess.run(
        [sys.executable, "scripts/reproduce_all.py", "--mode", "analysis", "--allow-hold", "--fail-on-missing"],
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 2


def test_reproduce_all_entrypoints_match_producer_registry_exactly():
    with (ROOT / "manifests" / "producer_status.tsv").open(encoding="utf-8") as handle:
        registered = sorted(row["entrypoint"].split("/")[-1] for row in csv.DictReader(handle, delimiter="\t"))
    assert sorted(PRODUCERS) == registered


def test_reproduce_all_default_surfaces_expected_value_mismatch(monkeypatch, tmp_path):
    monkeypatch.setattr(reproduce_all, "ROOT", tmp_path)
    monkeypatch.setattr(reproduce_all, "run_shared_chains", lambda **_: ({}, []))
    monkeypatch.setattr(reproduce_all, "all_object_ids", lambda _: ["Figure3"])

    def _mismatch_adapter(object_id, *, shared_dir, output_dir, repo_root):
        output_dir.mkdir(parents=True, exist_ok=True)
        (output_dir / f"{object_id}.status.json").write_text(json.dumps({"status": "HOLD_EXPECTED_VALUE_MISMATCH"}) + "\n", encoding="utf-8")
        return 2

    monkeypatch.setattr(reproduce_all, "run_object_adapter", _mismatch_adapter)
    assert reproduce_all.main(["--mode", "analysis", "--allow-hold"]) == 1


def test_reproduce_all_default_surfaces_malformed_shared_output(monkeypatch, tmp_path):
    monkeypatch.setattr(reproduce_all, "ROOT", tmp_path)
    monkeypatch.setattr(reproduce_all, "run_shared_chains", lambda **_: ({}, []))
    monkeypatch.setattr(reproduce_all, "all_object_ids", lambda _: ["Figure3"])

    def _malformed_adapter(object_id, *, shared_dir, output_dir, repo_root):
        shared_dir.mkdir(parents=True, exist_ok=True)
        (shared_dir / "B_PUBLIC_NUMERIC_OUTPUT.json").write_text("{not-json\n", encoding="utf-8")
        output_dir.mkdir(parents=True, exist_ok=True)
        (output_dir / f"{object_id}.status.json").write_text(
            json.dumps({"status": "HOLD_COMPONENT_NOT_AVAILABLE", "shared_chain": "B"}) + "\n",
            encoding="utf-8",
        )
        return 2

    monkeypatch.setattr(reproduce_all, "run_object_adapter", _malformed_adapter)
    assert reproduce_all.main(["--mode", "analysis", "--allow-hold"]) == 1
