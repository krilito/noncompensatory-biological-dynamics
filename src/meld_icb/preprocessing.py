"""Configuration-driven preprocessing with fail-closed coverage."""

from __future__ import annotations

from collections import defaultdict
from typing import Iterable, Mapping, Sequence

from .axes import available_gene_mean, sample_z


def collapse_duplicate_symbols(records: Iterable[tuple[str, float]]) -> dict[str, float]:
    """Collapse repeated gene symbols by arithmetic mean before axis scoring."""
    values: dict[str, list[float]] = defaultdict(list)
    for symbol, value in records:
        values[str(symbol)].append(float(value))
    if not values:
        raise ValueError("no gene records")
    return {symbol: sum(items) / len(items) for symbol, items in values.items()}


def axis_matrix(
    samples: Sequence[Mapping[str, float]],
    axis_genes: Mapping[str, Iterable[str]],
    *,
    ddof: int,
    minimum_coverage: float = 0.5,
) -> dict[str, list[float]]:
    if not samples:
        raise ValueError("at least one sample is required")
    raw: dict[str, list[float]] = {}
    for axis, genes in axis_genes.items():
        declared = list(genes)
        if not declared:
            raise ValueError(f"axis {axis} has no declared genes")
        rows: list[float] = []
        for sample in samples:
            observed = sum(1 for gene in declared if gene in sample and sample[gene] is not None)
            if observed / len(declared) < minimum_coverage:
                raise ValueError(f"axis {axis} fails minimum coverage")
            rows.append(available_gene_mean(sample, declared))
        raw[axis] = rows
    return {axis: sample_z(values, ddof=ddof) for axis, values in raw.items()}
