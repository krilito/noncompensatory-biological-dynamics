# PAIR Longitudinal Dataset

Unified, versioned longitudinal cancer dataset derived from the PAIR corpus:
patient → samples → longitudinal interval → treatment → clinical endpoint → expression references → provenance.

Main entry point: `releases/v0.1/PAIR_LONGITUDINAL_MASTER.parquet` (also `.csv.gz`).

## Build

```bash
python src/build_pair_dataset.py            # writes releases/v0.1/
```

The build is deterministic and read-only with respect to the corpus: it projects
`../corpus` canonical parquet tables plus the recovery pipeline in
`../corpus/src/patient_archive.py` (GIDE author-identity repair, CRS/lesion/RECIST-binary
recovery, native-semantics resolution). No second truth source is created.

## Layout

```
src/build_pair_dataset.py    # master-table builder (row = patient × interval × treatment)
src/normalization.py         # phase / treatment / endpoint / harmonization mappings
releases/v0.1/               # generated release (committed for review)
README.md  DATASET_CARD.md  SCHEMA.md
```

## Key invariants

- Native clinical labels are never rewritten; harmonization is an additive field.
- Patient identity merges only on confirmed evidence (author subject IDs,
  `global_patient_identity` confirmed roots). GSE20181/GSE5462 is a same-study
  resource duplicate whose patients are NOT merged.
- No globally batch-corrected expression matrix is produced; expression files are
  referenced, not embedded.
- NeoTRIP is present only as an access-status record (see LABEL_RECOVERY_LEDGER.csv).

See `DATASET_CARD.md` for scope/limitations and `SCHEMA.md` for column definitions.

## Repository position

In this public repository the `dataset/` directory publishes the generated v0.1
release tables, the deterministic builder (`src/`) and documentation for review.
Regenerating the release additionally requires the internal PAIR longitudinal
corpus (canonical parquet tables + patient-archive recovery pipeline); see
`docs/DATA_ACCESS.md` for how the underlying public resources are acquired.
