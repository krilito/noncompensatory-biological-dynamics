"""Axis arithmetic for available-gene means and sample standardisation."""

from __future__ import annotations

from math import sqrt
from typing import Iterable, Mapping


def available_gene_mean(values: Mapping[str, float], genes: Iterable[str]) -> float:
    observed = [float(values[g]) for g in genes if g in values and values[g] is not None]
    if not observed:
        raise ValueError("axis has no available genes")
    return sum(observed) / len(observed)


def sample_z(values: Iterable[float], ddof: int = 1) -> list[float]:
    xs = [float(v) for v in values]
    if len(xs) <= ddof:
        raise ValueError("insufficient values for requested ddof")
    mean = sum(xs) / len(xs)
    variance = sum((x - mean) ** 2 for x in xs) / (len(xs) - ddof)
    if variance == 0:
        raise ValueError("zero-variance axis is not scorable")
    scale = sqrt(variance)
    return [(x - mean) / scale for x in xs]
