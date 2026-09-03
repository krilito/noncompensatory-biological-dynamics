# Paired measurements and biological change

What paired biological measurements identify about biological change.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21975986.svg)](https://doi.org/10.5281/zenodo.21975986)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.11](https://img.shields.io/badge/python-3.11-blue.svg)](pyproject.toml)

![Melanoma tumour–immune microenvironment](docs/assets/readme_melanoma_tme.gif)

Paired measurements reveal biological change, but second-state content, recorded orientation, movement, task relevance, biological carrier, and interpretation domain are distinct inferential objects.

Melanoma immune checkpoint blockade. Paired pretreatment and on-treatment molecular measurements. When a biological state moves, what does that observed change actually identify?

## The question

A pretreatment and an on-treatment biopsy make change visible. But what exactly has been identified?

These interfaces are separable. They are not a composite score.

- second-state information
- recorded orientation
- coherent movement
- response specificity / incremental value
- biological carrier
- interpretation domain

## What we test

```text
Paired biological measurements
            ↓
     shared T/E/X/C state
            ↓
      PRE → ON movement
            ↓
What does the observed change identify?
```

![Paired PRE and ON states in T/E/X/C space](docs/assets/readme_texc_displacement.gif)

A/U/O then separates pretreatment state from pair content and recorded orientation, holding patients, folds and the learner fixed:

- **A** — pretreatment state
- **U** — complete unordered pair content
- **O** — the same pair with recorded orientation restored

$$d = a_{\mathrm{on}}-a_{\mathrm{pre}},\qquad m=(a_{\mathrm{pre}}+a_{\mathrm{on}})/2$$

$$A=a_{\mathrm{pre}},\qquad U=[m,\operatorname{vech}(dd^{\mathsf T})],\qquad O=[U,d]$$

MELD-ICB is the frozen four-axis evaluation object, not a clinical predictor or a complete biological world model.

## Main findings

### Pair content and orientation

In 27 paired patients, mean held-out loss was \(L_A=0.728\), \(L_U=0.640\), \(L_O=0.533\). The orientation-associated gain was \(G_{OU}=0.107\), with a grouped-patient bootstrap 95% interval of \(-0.354\) to \(0.729\).

Orientation-sensitive predictive structure was observed. Its population magnitude remained imprecise.

### Biological state moved

In PRJEB23709, 16 therapy records from 15 patients moved strongly and coherently: boundary score increased in 14/16 (paired Wilcoxon \(P=0.0017\)); cohort-mean cosine \(0.987\); median record-level cosine \(0.699\).

The sampled transcriptomic state moved across the observed treatment interval. Attribution specifically to ICB requires an untreated longitudinal comparator, which was not available.

### Movement was not equivalent to response

Boundary-change AUC was \(0.327\); transition-cosine AUC was \(0.364\). Adding observed change increased mean held-out loss (mean \(\Delta L=-0.0813\)).

The dominant displacement extended across response groups and did not improve mean held-out prediction beyond pretreatment state.

### Coordinates are not biological carriers

In CD45+-restricted material the tumour-proliferation coordinate remained numerically evaluable while \(86.293\%\) of its signal was lymphoid-carried and \(0.079\%\) malignant-carried (three cells).

A mathematical coordinate can remain computable after its intended biological carrier has changed.

## Context transfer

GSE78220 pretreatment melanoma AUC \(0.476\). IMvigor210 AUC \(0.400\). Selected axis-level phenotype structure could persist even when original response orientation did not.

Representation structure can persist while response interpretation fails to transfer.

> Paired observations identify change upstream of clinical benefit, biological referent and portable interpretation. Each stronger claim requires evidence of its own.

Statuses are non-compensatory. Full ceiling: [`docs/SCIENTIFIC_CLAIM_BOUNDARY.md`](docs/SCIENTIFIC_CLAIM_BOUNDARY.md).

## Observation interfaces

Frozen separator: \(-0.676T+0.173E+0X+1.135C-0.666\), threshold \(0.735\). Higher scores remain responder-like. AUCs below \(0.5\) are not inverted.

| Interface | Status |
|-----------|--------|
| B0 measurement integrity | Partial / historical preprocessing lineage **conflict** |
| B1 core state transfer (PRJEB23709, MORRISON-1-public) | **Supported** |
| B1 MGH biopsy-level transfer | **Not supported** |
| B2 movement | **Supported** |
| B3 outcome specificity | **Not supported** |
| B4 incremental information | **Not supported** |
| B5 biological carrier | Descriptive |
| B6 interpretation domain | Envelope declared |
| Intervention specificity | **Not tested** |

## Repository

| Path | Role |
|------|------|
| [`src/`](src/) | Public analysis implementation |
| [`configs/`](configs/) | Frozen contracts and cohort configs |
| [`scripts/`](scripts/) | Public entrypoints |
| [`source_data/`](source_data/) | Derived numeric Source Data |
| [`manifests/`](manifests/) | Dataset, producer, and numeric-reconciliation contracts |
| [`tests/`](tests/) | Unit, invariant, leakage, and reproduction tests |
| [`docs/`](docs/) | Claim boundary, data access, reproducibility scope |
| [`figures/reproduced/`](figures/reproduced/) | Quantitative redraws from Source Data |
| [`data/`](data/) | Local accession layout only; raw matrices are not shipped |

Figure 1 is conceptual author artwork and is not in this tree. Quantitative Figures 2–5 redraw from `source_data/`; publication Adobe boards are outside the code-reproduction claim.

## Reproduce

Python 3.11. Dependencies are lower-bounded in [`pyproject.toml`](pyproject.toml). There is no lockfile.

```text
python -m venv .venv
```

Windows: `.\.venv\Scripts\Activate.ps1`  
Linux / macOS: `source .venv/bin/activate`

```text
python -m pip install -e ".[dev]"
python -m pytest -q
python scripts/scan_public_repository.py
```

That is LEVEL 1: installation, frozen constants, invariants, and HOLD behavior. It does not recompute the paper from raw transcriptomes. GitHub Actions runs LEVEL 1 plus a Source Data redraw; a green badge is not full scientific reproduction.

LEVEL 2, from committed Source Data:

```text
python scripts/reproduce_core_results.py
python scripts/plot_figures.py
```

LEVEL 3 rebuilds from public accessions after `python scripts/acquire_data.py`. Missing inputs exit HOLD. They are never replaced by copied publication values. See [`docs/REPRODUCIBILITY_SCOPE.md`](docs/REPRODUCIBILITY_SCOPE.md).

Conda alternative: `conda env create -f environment.yml`.

## Data

No new raw sequencing was generated. Study inputs are public transcriptomic accessions. Committed `source_data/` holds derived numeric tables only. Raw matrices, single-cell counts, and restricted clinical objects are not redistributed.

See [`docs/DATA_ACCESS.md`](docs/DATA_ACCESS.md) and [`manifests/datasets.tsv`](manifests/datasets.tsv).

| Dataset | Role | In this repo |
|---------|------|----------------|
| Derived Source Data CSVs | Figure and table numbers | Derived tables |
| GSE91061, PRJEB23709, MORRISON-1-public, MGH | Analysis inputs | Accession only |
| GSE120575 (Sade-Feldman) | Carrier grounding | Accession only |
| GSE78220 | Envelope primary | Accession only |
| IMvigor210 | Comparison only | Accession only |

A public accession is not a redistribution license.

## Citation

This is a software/archive DOI, not a journal article DOI. The manuscript title is *What paired biological measurements identify about biological change*. Use [`CITATION.cff`](CITATION.cff). Software author: Yizhang Yang. Paper: Yizhang Yang (first); Chao Yang (corresponding, yangchao@qdu.edu.cn).

Release `v1.0.0-preprint` is archived at [10.5281/zenodo.21975986](https://doi.org/10.5281/zenodo.21975986). Concept DOI: [10.5281/zenodo.21975985](https://doi.org/10.5281/zenodo.21975985). No journal citation is invented here.

```bibtex
@software{yang2026noncompensatory,
  author = {Yang, Yizhang},
  title = {Non-compensatory Biological Dynamics},
  year = {2026},
  version = {v1.0.0-preprint},
  doi = {10.5281/zenodo.21975986},
  url = {https://github.com/krilito/noncompensatory-biological-dynamics},
  note = {Accompanying manuscript: What paired biological measurements identify about biological change. Corresponding author: Chao Yang.}
}
```

## License

Code is [MIT](LICENSE). Source datasets remain under their original terms; MIT does not relicense them.
