"""Small dependency-free statistics used by public producers."""

from __future__ import annotations

from typing import Iterable, Sequence


def auc(scores: Sequence[float], labels: Sequence[int]) -> float:
    if len(scores) != len(labels) or not scores:
        raise ValueError("scores and labels must have equal nonzero length")
    positives = [s for s, y in zip(scores, labels) if int(y) == 1]
    negatives = [s for s, y in zip(scores, labels) if int(y) == 0]
    if not positives or not negatives:
        raise ValueError("AUC requires both response classes")
    wins = sum(1.0 if p > n else 0.5 if p == n else 0.0 for p in positives for n in negatives)
    return wins / (len(positives) * len(negatives))


def binary_log_loss(probabilities: Iterable[float], labels: Iterable[int]) -> float:
    import math

    ps = [min(max(float(p), 1e-15), 1 - 1e-15) for p in probabilities]
    ys = [int(y) for y in labels]
    if len(ps) != len(ys) or not ps:
        raise ValueError("probabilities and labels must have equal nonzero length")
    return -sum(y * math.log(p) + (1 - y) * math.log(1 - p) for p, y in zip(ps, ys)) / len(ps)


def registered_auc_inference(
    scores: Sequence[float],
    labels: Sequence[int],
    *,
    identifiers: Sequence[str],
    identifier_order: str,
    seed: int,
    bootstrap_resamples: int,
    permutation_resamples: int,
    inference_contract: str,
    rng_algorithm: str,
    stream_policy: str,
    permutation_comparator: str,
    correction: str,
) -> dict[str, object]:
    """Execute one of the two closed A/D inference contracts.

    Each contract has one exact resample/seed/RNG/stream combination and every
    mismatch fails closed. Stream sharing and permutation-buffer persistence
    are derived from the contract; they are not public switches.
    """
    import random

    if len(scores) != len(labels) or len(scores) != len(identifiers):
        raise ValueError("registered inference requires aligned scores, labels, and identifiers")
    if not scores or len(set(int(label) for label in labels)) < 2:
        raise ValueError("registered inference requires aligned two-class observations")
    contract = str(inference_contract).strip()
    algorithm = str(rng_algorithm).strip()
    expected = {
        "A_EXT_AUC_PERMUTATION_V1": {
            "bootstrap_resamples": 2000,
            "permutation_resamples": 1000,
            "seed": 42,
            "rng_algorithm": "python_random_v1",
            "shared_rng_stream": False,
            "persistent_permutation_buffer": False,
            "stream_policy": "independent_bootstrap_and_permutation",
        },
        "D_OPERATING_ENVELOPE_REGISTERED_V1": {
            "bootstrap_resamples": 10000,
            "permutation_resamples": 10000,
            "seed": 20260731,
            "rng_algorithm": "python_random_shared_v1",
            "shared_rng_stream": True,
            "persistent_permutation_buffer": True,
            "stream_policy": "shared_bootstrap_then_permutation_persistent_buffer",
        },
    }.get(contract)
    if expected is None:
        raise ValueError(f"unsupported inference contract: {contract}")
    declared_stream_policy = str(stream_policy).strip()
    observed = {
        "bootstrap_resamples": int(bootstrap_resamples),
        "permutation_resamples": int(permutation_resamples),
        "seed": int(seed),
        "rng_algorithm": algorithm,
        "shared_rng_stream": expected["shared_rng_stream"],
        "persistent_permutation_buffer": expected["persistent_permutation_buffer"],
        "stream_policy": declared_stream_policy,
    }
    for key, value in expected.items():
        if observed[key] != value:
            raise ValueError(f"{contract} requires {key}={value!r}, observed {observed[key]!r}")
    if permutation_comparator != "inclusive_absolute_deviation_ge":
        raise ValueError("only inclusive_absolute_deviation_ge is legal")
    if correction != "plus_one":
        raise ValueError("only plus_one correction is legal")
    normalized_ids = [str(identifier).strip() for identifier in identifiers]
    if any(not identifier for identifier in normalized_ids) or len(set(normalized_ids)) != len(normalized_ids):
        raise ValueError("registered inference requires nonempty unique identifiers")

    def patient_key(identifier: str) -> tuple[int, int | str, str]:
        try:
            parsed = int(identifier)
        except ValueError:
            return (1, identifier, identifier)
        if str(parsed) == identifier:
            return (0, parsed, identifier)
        return (1, identifier, identifier)

    if identifier_order == "sample_id":
        order = sorted(range(len(normalized_ids)), key=lambda index: normalized_ids[index])
    elif identifier_order == "patient_id":
        order = sorted(range(len(normalized_ids)), key=lambda index: patient_key(normalized_ids[index]))
    else:
        raise ValueError("identifier_order must be sample_id or patient_id")
    ordered_scores = [float(scores[index]) for index in order]
    ordered_labels = [int(labels[index]) for index in order]
    ordered_ids = [normalized_ids[index] for index in order]
    positives = [index for index, label in enumerate(ordered_labels) if label == 1]
    negatives = [index for index, label in enumerate(ordered_labels) if label == 0]
    bootstrap_rng = random.Random(seed)
    permutation_rng = bootstrap_rng if expected["shared_rng_stream"] else random.Random(seed)
    bootstrap: list[float] = []
    for _ in range(bootstrap_resamples):
        draw = [bootstrap_rng.choice(positives) for _ in positives] + [bootstrap_rng.choice(negatives) for _ in negatives]
        bootstrap.append(auc([ordered_scores[index] for index in draw], [ordered_labels[index] for index in draw]))
    ordered = sorted(bootstrap)
    low = ordered[max(0, int(0.025 * len(ordered)) - 1)]
    high = ordered[min(len(ordered) - 1, int(0.975 * len(ordered)))]
    observed = auc(ordered_scores, ordered_labels)
    exceed = 0
    permutation_buffer = list(ordered_labels)
    for _ in range(permutation_resamples):
        # A copies canonical labels for every draw.  D intentionally mutates
        # one persistent buffer under its separate registered contract.
        permuted = permutation_buffer if expected["persistent_permutation_buffer"] else list(ordered_labels)
        permutation_rng.shuffle(permuted)
        exceed += int(abs(auc(ordered_scores, permuted) - 0.5) >= abs(observed - 0.5))
    return {
        "auc": observed,
        "ci_95": (low, high),
        "p_two_sided": (exceed + 1) / (permutation_resamples + 1),
        "bootstrap_resamples": bootstrap_resamples,
        "permutation_resamples": permutation_resamples,
        "seed": seed,
        "rng_algorithm": algorithm,
        "stream_policy": expected["stream_policy"],
        "identifier_order": identifier_order,
        "ordered_identifiers": tuple(ordered_ids),
        "permutation_comparator": permutation_comparator,
        "correction": correction,
        "inference_contract": contract,
    }
