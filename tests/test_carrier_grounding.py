from pathlib import Path
import os

import pandas as pd
import pytest

from meld_icb.carrier_grounding import CarrierInputError, _assign_carrier_lineage, run_carrier_grounding


ROOT = Path(__file__).resolve().parents[1]
_TRUTH_ENV = os.environ.get("MELD_TRUTH_ROOT")
TRUTH = Path(_TRUTH_ENV) if _TRUTH_ENV else Path("__missing_truth__")


def test_zero_signal_cd8_cd4_tie_uses_source_cd4_fallback():
    expr = pd.DataFrame({
        "c_zero": [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0],
    }, index=["CD3D", "CD3E", "CD3G", "CD8A", "CD8B", "CD4", "IL7R"])
    assert _assign_carrier_lineage(expr)["c_zero"] == "CD4_T"


def test_nonzero_cd8_cd4_tie_is_fail_closed():
    expr = pd.DataFrame({
        "c_tie": [1.0, 1.0, 1.0, 0.5, 0.5, 0.5, 0.5],
    }, index=["CD3D", "CD3E", "CD3G", "CD8A", "CD8B", "CD4", "IL7R"])
    with pytest.raises(CarrierInputError, match="ambiguity"):
        _assign_carrier_lineage(expr)


def test_missing_carrier_matrix_is_hold(tmp_path):
    result = run_carrier_grounding(truth_root=tmp_path / "missing", config_path=ROOT / "configs" / "carrier_grounding.yaml")
    assert result["verdict"] == "HOLD"


@pytest.mark.skipif(not TRUTH.exists(), reason="truth repository is not available in a fresh clone")
def test_truth_carrier_gate_and_source_summary():
    result = run_carrier_grounding(truth_root=TRUTH, config_path=ROOT / "configs" / "carrier_grounding.yaml")
    assert result["verdict"] == "GATES_PASSED"
    assert result["n_cells"] == 10363
    assert result["gates"]["A_lymphoid_share"]["observed"] == 0.8629283905029297
    t_axis = next(row for row in result["summary"] if row["signature"] == "T_axis")
    assert t_axis["lymphoid_share"] == pytest.approx(0.8629283383488655, abs=1e-15)
    assert t_axis["malignant_share"] == pytest.approx(0.0007919137715362012, abs=1e-18)
    assert result["semantic_limits"]["patient_ci"] == "NOT_ESTIMABLE"
