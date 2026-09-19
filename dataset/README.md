# PAIR Longitudinal Dataset

Unified, versioned longitudinal cancer dataset derived from the PAIR corpus:
patient → samples → longitudinal interval → treatment → clinical endpoint → expression references → provenance.

Main entry point: `releases/v0.1.2/PAIR_LONGITUDINAL_MASTER.parquet` (also `.csv.gz`).

Reading the release tables requires `pyarrow` (`pip install "pyarrow>=14"`, or
`pip install -e ".[dev]"` from the repository root); every published table is parquet, and
`dataset/tests` reads them. Path locators in the release are project-relative POSIX, never
machine-absolute - see the portability section of `SCHEMA.md`.

## Build

```bash
python src/build_pair_dataset.py            # writes releases/v0.1.2/
python src/build_gene_space.py              # writes gene_space/ (v0.2.1 canonical gene space)
python src/build_expression_layer.py        # writes expression_layer/ + local C1 matrices
```

The build is deterministic and read-only with respect to the corpus: it projects
`../corpus` canonical parquet tables plus the recovery pipeline in
`../corpus/src/patient_archive.py` (GIDE author-identity repair, CRS/lesion/RECIST-binary
recovery, native-semantics resolution). No second truth source is created.

## Layout

```
src/build_pair_dataset.py    # master-table builder (row = patient × interval × treatment)
src/normalization.py         # phase / treatment / endpoint / harmonization mappings
src/build_gene_space.py       # level B: HGNC canonical gene space + feature maps
src/expression_transforms.py  # level C1: the declared-scale transformation contract
src/build_expression_layer.py # level C1: per-source quantitative matrices + sample binding
releases/v0.1.2/             # generated release (committed for review)
gene_space/                  # level B artifacts (published)
expression_layer/            # level C1 contracts, metrics, manifests (published)
README.md  DATASET_CARD.md  SCHEMA.md
```

## Key invariants

- Native clinical labels are never rewritten; harmonization is an additive field.
- Patient identity merges only on confirmed evidence (author subject IDs,
  `global_patient_identity` confirmed roots, shared GSM sample accessions with
  agreeing native patient ids). GSE20181/GSE5462 patients are merged at
  `patient_uid` level (v0.1.2) to prevent cross-resource train/test leakage; both
  accession-specific resource records are retained.
- No globally batch-corrected expression matrix is produced; expression files are
  referenced, not embedded. The C1 matrices in `expression_layer/` are per-source and
  per-scale only: nothing is z-scored, quantile-normalized jointly, rank-transformed or
  batch-corrected across cohorts, and no derived matrix is committed.
- File locators in published artifacts are project-relative POSIX, never machine
  absolute; `tests/test_master_invariants.py` enforces this over every release and layer
  table. See the "Path locators and portability" section of `SCHEMA.md`.
- NeoTRIP is present only as an access-status record (see LABEL_RECOVERY_LEDGER.csv).

See `DATASET_CARD.md` for scope/limitations and `SCHEMA.md` for column definitions.
