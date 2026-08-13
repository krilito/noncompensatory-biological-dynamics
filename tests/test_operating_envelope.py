import pytest

from meld_icb.operating_envelope import classify_auc_interval


def test_declared_near_null_and_unresolved_are_not_inferred_from_ci():
    assert classify_auc_interval(0.571, 0.371, 0.764, declared_status="NEAR_NULL") == "NEAR_NULL"
    assert classify_auc_interval(0.411, 0.173, 0.667, declared_status="UNRESOLVED") == "NEAR_NULL"
    assert classify_auc_interval(0.400, 0.327, 0.476, declared_status="INVERTED") == "INVERTED"
    assert classify_auc_interval(None, None, None, declared_status="HOLD", hold=True) == "HOLD"
    assert classify_auc_interval(0.411, 0.173, 0.667, declared_status="INVERTED") == "NEAR_NULL"
