from pathlib import Path

from meld_icb.config import load_freeze
from scripts.reproduce_core_results import load_contract, reconcile


ROOT = Path(__file__).resolve().parents[1]


def test_source_data_matches_manuscript_contract():
    rows = reconcile(ROOT)
    blockers = [row["RESULT"] for row in rows if row["STATUS"] != "PASS"]
    assert rows, "numeric contract is empty"
    assert not blockers, blockers


def test_frozen_boundary_matches_contract_and_machine_config():
    contract = load_contract(ROOT)["frozen_boundary"]
    freeze = load_freeze(ROOT)
    assert freeze["boundary"]["coefficients"] == {
        "T": contract["T"],
        "E": contract["E"],
        "X": contract["X"],
        "C": contract["C"],
    }
    assert freeze["boundary"]["intercept"] == contract["intercept"]
    assert freeze["boundary"]["threshold"] == contract["threshold"]


def test_historical_lineage_columns_are_retained():
    contract = load_contract(ROOT)
    header = (ROOT / "source_data" / "figure2_auc.csv").read_text(encoding="utf-8").splitlines()[0]
    for column in contract["historical_lineage"]["retained_columns"]:
        assert column in header
