# Non-compensatory evaluation of learned biological dynamics

This repository accompanies the preprint

> **Non-compensatory evaluation of state, movement and meaning in learned
> biological dynamics**

MELD-ICB is used as a real-data evaluation object. The repository does not
claim a complete biological world model, a clinical predictor, or a causal
intervention effect.

Software and code: Yizhang Yang.

Paper: Yizhang Yang (first author); Chao Yang (corresponding author,
yangchao@qdu.edu.cn).

## Evidence status

Statuses are non-compensatory. A supported result does not rescue an
unsupported, descriptive, conflict, or not-tested result.

| Component | Status |
|-----------|--------|
| B0 measurement integrity | Partial / historical preprocessing lineage **conflict** |
| B1 core state transfer (PRJEB23709, MORRISON-1-public) | **Supported** |
| B1 MGH biopsy-level transfer | **Not supported** |
| B1 MGH patient-level sensitivity | Supportive (different analysis unit) |
| B2 movement existence | **Supported** |
| B3 outcome specificity | **Not supported** |
| B4 incremental information | **Not supported** |
| B5 semantic grounding | Descriptive |
| B6 operating envelope | Envelope declared |
| Intervention specificity | **Not tested** |

The public code has one canonical preprocessing path,
`PRE_MEAN_THEN_AXIS_Z` with `ddof=1`. The B0 conflict is a historical
lineage fact retained in `source_data/figure2_auc.csv`
(`auc_canonical` vs `auc_historical`). It is not a second live
preprocessor that a user can switch on.

Frozen separator: `-0.676 T + 0.173 E + 0 X + 1.135 C - 0.666`,
threshold `0.735`. Higher scores remain responder-like. AUCs below 0.5
are not inverted.

## Install

Python 3.11 is the registered runtime.

```text
git clone https://github.com/krilito/noncompensatory-biological-dynamics.git
cd noncompensatory-biological-dynamics
python -m venv .venv
```

Windows PowerShell:

```text
.\.venv\Scripts\Activate.ps1
python -m pip install -U pip
python -m pip install -e ".[dev]"
```

Linux / macOS:

```text
source .venv/bin/activate
python -m pip install -U pip
python -m pip install -e ".[dev]"
```

Conda alternative: `conda env create -f environment.yml`.

Dependencies are lower-bounded in `pyproject.toml`. There is no lockfile
and no private or editable local dependency.

## Reproduction levels

These levels are not interchangeable.

### LEVEL 1 — tests only

Verifies installation, frozen constants, scientific invariants, repository
boundary, and missing-data HOLD behavior. It does **not** recompute the
paper from raw transcriptomes.

```text
python -m pytest -q
python scripts/scan_public_repository.py
```

GitHub Actions runs this level plus a Source Data figure redraw. The CI
badge, if added, means only that. It does not mean full scientific
reproduction.

### LEVEL 2 — reproduce from committed Source Data

Reconciles preprint numbers with committed tables and redraws quantitative
Figures 2–5.

```text
python scripts/reproduce_core_results.py
python scripts/plot_figure2.py
python scripts/plot_figure3.py
python scripts/plot_figures.py
```

This redraws numbers and quantitative panels. It is **not** a rebuild of
raw GEO/ENA matrices, and it is **not** a pixel copy of Adobe publication
boards.

### LEVEL 3 — rebuild from public accession data

Requires the user to obtain source datasets under their original terms
and point a truth root at the expected layout.

```text
python scripts/acquire_data.py
python scripts/recompute_results.py --truth-root <path>
python scripts/reproduce_all.py --mode analysis --truth-root <path>
```

Missing inputs exit with an explicit HOLD. They are never replaced by
copied publication values.

## Figures

| Figure | Public reproduction |
|--------|---------------------|
| Fig. 1 | Conceptual author artwork. Not code-generated. |
| Fig. 2 | Quantitative panels from `source_data/figure2*.csv`. |
| Fig. 3 | Quantitative panels from `source_data/figure3*.csv`. |
| Fig. 4 | Quantitative panels from `source_data/figure4*.csv`. Publication board is an Adobe master. |
| Fig. 5 | Panels a, b and d from committed `source_data/figure5*.csv` summaries. Panel c is accession-only: IMvigor210 sample-level scores are not redistributed, so LEVEL 2 does not redraw that strip plot. Publication board is an Adobe master. |

Color is not a scientific authority. Numeric values, sample sizes, AUCs,
intervals, and status labels are.

## Data

Committed `source_data/` tables are derived numeric Source Data classified
for redistribution. Raw expression matrices, single-cell count matrices,
restricted clinical objects, and private identifier maps are **not**
included and must not be added.

A public accession is not a redistribution license. See
`docs/DATA_ACCESS.md` and `manifests/datasets.tsv`.

| Dataset | Role | Redistribution |
|---------|------|----------------|
| Derived Source Data CSVs | Figure/table numbers | Derived tables only |
| GSE91061, PRJEB23709, MORRISON-1-public, MGH GEO studies | Analysis inputs | Accession only |
| GSE120575 (Sade-Feldman) | Carrier grounding | Accession only |
| GSE78220 | Envelope primary | Accession only |
| IMvigor210 | Comparison only | Accession only |
| TCGA-SKCM, GSE273583, GSE123813 | Context | Accession only |
| GSE282471 | Superseded hold | Accession only; not an active claim |

## Repository map

- `src/` public analysis implementation
- `configs/` frozen contracts
- `scripts/` public entrypoints
- `source_data/` derived numeric tables
- `manifests/` dataset, producer, and numeric-reconciliation contracts
- `tests/` unit, invariant, leakage, and reproduction tests
- `docs/` claim boundary, data access, and reproducibility scope

## Citation, license, version

Use `CITATION.cff`. The software author is Yizhang Yang. The paper
citation lists both authors and names Chao Yang as corresponding author.
Code is MIT. Source datasets remain under their original terms; MIT does
not relicense them.

Intended public tag for the arXiv v1 companion is `v1.0.0-preprint`.
The arXiv identifier and Zenodo DOI are placeholders until those records
exist. Do not invent them.

## Limitations

- B3 is a small-sample null (5 non-responder records), not proof that
  movement is generally response-independent.
- B4 mean ΔL is negative; adding change did not improve held-out loss.
- B5 is pooled-cell and descriptive. It is not a matched bulk comparator
  and not a causal substrate switch.
- Intervention specificity was not tested.
- Cohort-relative normalization is part of the locked measurement system.
- Figure 1 and the Adobe masters for Figures 4 and 5 are outside the
  code-reproduction claim.

Full claim ceiling: `docs/SCIENTIFIC_CLAIM_BOUNDARY.md`.
