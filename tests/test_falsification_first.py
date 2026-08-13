from __future__ import annotations

import math

from meld_icb.carrier_grounding import established_signature_gene_sets
from meld_icb.falsification_statistics import (
    benjamini_hochberg,
    change_coupling,
    exact_auc_permutation_p,
    paired_change_control,
    state_movement_metrics,
)
from meld_icb.falsification_first import _overlap_audit


def test_established_signatures_are_project_locked() -> None:
    signatures = established_signature_gene_sets()
    assert signatures["CYT"] == ("GZMA", "PRF1")
    assert len(signatures["IFNG6"]) == 6
    assert len(signatures["TIS18"]) == 18


def test_state_movement_metrics_keep_channels_aligned() -> None:
    pre = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]
    on = [0.1, 0.2, 0.3, 1.3, 1.4, 1.5]
    labels = [0, 0, 0, 1, 1, 1]
    patients = [str(index) for index in range(6)]
    result = state_movement_metrics(pre, on, labels, patients, bootstrap_repeats=200, seed=7)
    assert result["on"]["auc"] == 1.0
    assert result["movement"]["auc"] == 1.0
    assert result["n_R"] == result["n_NR"] == 3
    assert result["bootstrap_repeats_valid"] > 0


def test_exact_permutation_is_two_sided() -> None:
    labels = [0, 0, 1, 1]
    high = exact_auc_permutation_p([0.0, 0.1, 0.9, 1.0], labels)
    low = exact_auc_permutation_p([1.0, 0.9, 0.1, 0.0], labels)
    assert high == low
    assert math.isclose(high, 2 / 6)


def test_positive_control_and_coupling_are_finite() -> None:
    pre = [0.0, 0.2, 0.4, 0.6]
    on = [0.8, 1.1, 1.5, 2.0]
    patients = ["a", "b", "c", "d"]
    control = paired_change_control(pre, on, patients, bootstrap_repeats=50, seed=3)
    coupling = change_coupling(pre, on)
    assert control["fraction_positive"] == 1.0
    assert control["median_change"] > 0
    assert math.isfinite(coupling["spearman_pre_vs_delta_rho"])


def test_benjamini_hochberg_is_monotone_in_rank() -> None:
    adjusted = benjamini_hochberg([0.01, 0.04, 0.03])
    assert adjusted[0] <= adjusted[2] <= adjusted[1]
    assert all(0 <= value <= 1 for value in adjusted)


def test_overlap_audit_does_not_treat_source_local_tokens_as_patient_matches(monkeypatch) -> None:
    monkeypatch.setattr(
        "meld_icb.falsification_first.read_gse91061_on_treatment_labels",
        lambda truth: [{"patient_id": "Pt1", "sample_id": "Pt1_On"}],
    )
    result = _overlap_audit(
        None,
        [{"patient_id": 1, "pre_sample_id": "PD1_1_PRE", "edt_sample_id": "PD1_1_EDT", "therapy_type": "PD1"}],
    )
    assert result["N_overlap_patients"] == 0
    assert result["unqualified_numeric_token_collisions"] == ["1"]
