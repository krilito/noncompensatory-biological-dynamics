# Reproducibility scope

The public package supports three distinct levels. They are not equivalent.

## LEVEL 1 — tests only

```text
python -m pytest -q
python scripts/scan_public_repository.py
```

This checks installation, frozen coefficients and threshold, scientific
invariants, repository-boundary rules, and missing-data HOLD behavior.
CI runs this level. It does not download public datasets and does not
claim full scientific reproduction.

## LEVEL 2 — committed Source Data

```text
python scripts/reproduce_core_results.py
python scripts/plot_figures.py
```

This reconciles preprint numbers with committed derived tables and redraws
quantitative Figures 2–5 except Figure 5 panel c. IMvigor210 sample-level
scores are not redistributed; panel c shows only the paper-reported
summary AUC/CI/P and the accession. LEVEL 2 does **not** rebuild
preprocessing from raw GEO/ENA/GDC matrices and does **not** reproduce
Adobe editorial boards pixel for pixel.

Figure 1 is conceptual author artwork and has no numeric producer.

## LEVEL 3 — accession rebuild

```text
python scripts/acquire_data.py
python scripts/recompute_results.py --truth-root <path>
python scripts/reproduce_all.py --mode analysis --truth-root <path>
```

This recomputes registered chains from separately acquired source
datasets. Missing data, schema, or provenance remain HOLD. Copied
publication values are not substituted.

The release is generated one way from a controlled upstream workbench.
Its Git history is independent and is not an input authority for frozen
values.
