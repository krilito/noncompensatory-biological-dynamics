import pytest

from meld_icb.carrier_grounding import cell_count_shares, weighted_signal_shares
from meld_icb.statistics import auc, registered_auc_inference


def test_auc_and_signal_share_are_computed_from_inputs():
    assert auc([0.1, 0.9, 0.2, 0.8], [0, 1, 0, 1]) == 1.0
    assert weighted_signal_shares({"lymphoid": 3.0, "myeloid": 1.0}) == {"lymphoid": 0.75, "myeloid": 0.25}
    assert cell_count_shares({"lymphoid": 3, "myeloid": 1}) == {"lymphoid": 0.75, "myeloid": 0.25}
    with pytest.raises(ValueError):
        weighted_signal_shares({"bad": -1.0})


def test_registered_inference_normalizes_identifier_order_and_streams():
    scores = [0.1, 0.9, 0.2, 0.8, 0.4, 0.7]
    labels = [0, 1, 0, 1, 0, 1]
    identifiers = ["s3", "s1", "s5", "s2", "s6", "s4"]
    baseline = registered_auc_inference(
        scores,
        labels,
        identifiers=identifiers,
        identifier_order="sample_id",
        seed=42,
        bootstrap_resamples=2000,
        permutation_resamples=1000,
        inference_contract="A_EXT_AUC_PERMUTATION_V1",
        rng_algorithm="python_random_v1",
        stream_policy="independent_bootstrap_and_permutation",
        permutation_comparator="inclusive_absolute_deviation_ge",
        correction="plus_one",
    )
    permutation = [4, 1, 5, 0, 3, 2]
    shuffled = registered_auc_inference(
        [scores[index] for index in permutation],
        [labels[index] for index in permutation],
        identifiers=[identifiers[index] for index in permutation],
        identifier_order="sample_id",
        seed=42,
        bootstrap_resamples=2000,
        permutation_resamples=1000,
        inference_contract="A_EXT_AUC_PERMUTATION_V1",
        rng_algorithm="python_random_v1",
        stream_policy="independent_bootstrap_and_permutation",
        permutation_comparator="inclusive_absolute_deviation_ge",
        correction="plus_one",
    )
    assert shuffled["auc"] == baseline["auc"]
    assert shuffled["ci_95"] == baseline["ci_95"]
    assert shuffled["p_two_sided"] == baseline["p_two_sided"]
    with pytest.raises(ValueError):
        registered_auc_inference(
            scores,
            labels,
            identifiers=identifiers,
            identifier_order="sample_id",
            seed=42,
            bootstrap_resamples=1999,
            permutation_resamples=1000,
            inference_contract="A_EXT_AUC_PERMUTATION_V1",
            rng_algorithm="python_random_v1",
            stream_policy="independent_bootstrap_and_permutation",
            permutation_comparator="inclusive_absolute_deviation_ge",
            correction="plus_one",
        )
    assert baseline["rng_algorithm"] == "python_random_v1"
    assert baseline["inference_contract"] == "A_EXT_AUC_PERMUTATION_V1"


def test_registered_contract_stream_policies_cannot_cross_use():
    scores = [0.1, 0.9, 0.2, 0.8, 0.4, 0.7]
    labels = [0, 1, 0, 1, 0, 1]
    identifiers = ["s1", "s2", "s3", "s4", "s5", "s6"]
    independent = registered_auc_inference(
        scores,
        labels,
        identifiers=identifiers,
        identifier_order="sample_id",
        seed=42,
        bootstrap_resamples=2000,
        permutation_resamples=1000,
        inference_contract="A_EXT_AUC_PERMUTATION_V1",
        rng_algorithm="python_random_v1",
        stream_policy="independent_bootstrap_and_permutation",
        permutation_comparator="inclusive_absolute_deviation_ge",
        correction="plus_one",
    )
    shared = registered_auc_inference(
        scores,
        labels,
        identifiers=identifiers,
        identifier_order="sample_id",
        seed=20260731,
        bootstrap_resamples=10000,
        permutation_resamples=10000,
        inference_contract="D_OPERATING_ENVELOPE_REGISTERED_V1",
        rng_algorithm="python_random_shared_v1",
        stream_policy="shared_bootstrap_then_permutation_persistent_buffer",
        permutation_comparator="inclusive_absolute_deviation_ge",
        correction="plus_one",
    )
    assert independent["inference_contract"] == "A_EXT_AUC_PERMUTATION_V1"
    assert independent["stream_policy"] == "independent_bootstrap_and_permutation"
    assert independent["rng_algorithm"] == "python_random_v1"
    assert shared["inference_contract"] == "D_OPERATING_ENVELOPE_REGISTERED_V1"
    assert shared["stream_policy"] == "shared_bootstrap_then_permutation_persistent_buffer"
    assert shared["rng_algorithm"] == "python_random_shared_v1"
    assert shared["p_two_sided"] != independent["p_two_sided"]


def test_registered_inference_rejects_declared_stream_hybrid():
    with pytest.raises(ValueError):
        registered_auc_inference(
            [0.1, 0.9, 0.2, 0.8],
            [0, 1, 0, 1],
            identifiers=["s1", "s2", "s3", "s4"],
            identifier_order="sample_id",
            seed=42,
            bootstrap_resamples=2000,
            permutation_resamples=1000,
            inference_contract="A_EXT_AUC_PERMUTATION_V1",
            rng_algorithm="python_random_v1",
            stream_policy="shared_bootstrap_then_permutation_persistent_buffer",
            permutation_comparator="inclusive_absolute_deviation_ge",
            correction="plus_one",
        )
