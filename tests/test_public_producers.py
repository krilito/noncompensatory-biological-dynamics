import json
from pathlib import Path

from meld_icb.producers import figure1, figure2, figure3, figure4, figure5
from scripts.public_runner import produce_computed


ROOT = Path(__file__).resolve().parents[1]


def test_retired_raw_producers_fail_closed_without_statistics():
    for producer in (figure1, figure2, figure3, figure4, figure5):
        assert producer("raw.tsv", ROOT, ROOT)["status"] == "HOLD_COMPONENT_NOT_AVAILABLE"


def test_public_runner_ignores_raw_inputs_and_projects_shared_chain_only(tmp_path):
    assert produce_computed("Figure1", ["C1"], ["GSE91061"], "figure1", ["raw.tsv"], output_dir=str(tmp_path)) == 2
    payload = json.loads((tmp_path / "Figure1.status.json").read_text(encoding="utf-8"))
    assert payload["status"] == "HOLD_RAW_RECOMPUTE_RETIRED"
    assert not (tmp_path / "Figure1.svg").exists()
