# Non-compensatory Biological Dynamics

Code and reproducibility materials for **“Non-compensatory evaluation of
state, movement and meaning in learned biological dynamics.”** The study asks
which evidence conditions must qualify before predicted biological change
warrants a learned-dynamics claim for future biological world models. Melanoma
immune checkpoint blockade is the empirical stress test.

## Reproducibility scope

This repository contains the Python implementation, frozen configurations,
tests, redistributable derived Source Data, acquisition contracts, and scripts
for the main quantitative results and quantitative figure redraws. Figure 1 is
conceptual author artwork. The scripts for Figures 2–5 reproduce quantitative
panels from committed tables; they do not reproduce the final Adobe editorial
composition pixel for pixel.

No raw expression matrix, restricted cohort object, private identifier map,
manuscript source, Adobe file, internal audit record, or internal process trace
is included. There is no decorative R port: the active producers are Python.

This repository is a one-way release artifact from a controlled research
workbench. Changes made here are not imported upstream and must not be used to
silently alter frozen scientific values or analysis contracts.

## Quick start

Python 3.11 is the registered runtime.

```text
python -m venv .venv
python -m pip install -e ".[dev]"
python -m pytest -q
python scripts/plot_figures.py
```

Quantitative redraws are written to `figures/reproduced/`. Run the complete
public workflow with:

```text
python scripts/reproduce_all.py --mode figures
```

Analysis recomputation requires source datasets obtained under their original
terms. See `docs/DATA_ACCESS.md`, `manifests/datasets.tsv`, and
`data/README.md`. Missing inputs produce explicit HOLD states; they are never
replaced by copied publication values or fabricated data.

## Repository map

- `src/`: analysis implementation.
- `configs/`: frozen model and cohort contracts.
- `scripts/`: public analysis and figure entrypoints.
- `source_data/`: redistributable derived numeric tables.
- `manifests/`: dataset, producer, and output provenance contracts.
- `tests/`: scientific regression, behavior, and release-boundary tests.
- `figures/`: generated-output destination and scope note.

## Citation and license

Use `CITATION.cff` for software citation metadata. Code is released under the
MIT License. Source datasets and derived tables remain subject to the original
data-source terms described in the manifests; the MIT License does not grant
rights to third-party datasets.
