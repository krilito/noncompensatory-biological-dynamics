from pathlib import Path

import pytest

from meld_icb.config import load_freeze
from meld_icb.preprocessing import axis_matrix, collapse_duplicate_symbols


ROOT = Path(__file__).resolve().parents[1]


def test_frozen_constants_are_machine_readable():
    freeze = load_freeze(ROOT)
    assert freeze["boundary"]["coefficients"] == {"T": -0.676, "E": 0.173, "X": 0.0, "C": 1.135}
    assert freeze["boundary"]["intercept"] == -0.666
    assert freeze["boundary"]["threshold"] == 0.735
    assert freeze["preprocessing"]["ddof"] == 1
    assert freeze["preprocessing"]["minimum_axis_coverage"] == 0.5


def test_axis_coverage_fail_closed_and_duplicate_mean_rule():
    samples = [{"A": 1.0, "A_duplicate": 3.0, "B": 2.0}, {"A": 2.0, "A_duplicate": 4.0, "B": 3.0}]
    collapsed = collapse_duplicate_symbols([("A", 1.0), ("A", 3.0), ("B", 2.0)])
    assert collapsed["A"] == 2.0
    values = axis_matrix(samples, {"T": ["A", "A_duplicate"], "E": ["B"]}, ddof=1, minimum_coverage=0.5)
    assert values["T"][0] < values["T"][1]
    with pytest.raises(ValueError, match="minimum coverage"):
        axis_matrix([{"A": 1.0}], {"T": ["A", "MISSING", "MISSING2"]}, ddof=1, minimum_coverage=0.5)
