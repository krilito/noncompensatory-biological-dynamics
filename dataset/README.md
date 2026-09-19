# PAIR Longitudinal Dataset

Unified, versioned longitudinal cancer dataset derived from the PAIR corpus:
patient → samples → longitudinal interval → treatment → clinical endpoint → expression references → provenance.

Main entry point: `releases/v0.1.2/PAIR_LONGITUDINAL_MASTER.parquet` (also `.csv.gz`).

## Build

```bash
python src/build_pair_dataset.py            # writes releases/v0.1.2/
```

The build is deterministic and read-only with respect to the corpus: it projects
`../corpus` canonical parquet tables plus the recovery pipeline in
`../corpus/src/patient_archive.py` (GIDE author-identity repair, CRS/lesion/RECIST-binary
recovery, native-semantics resolution). No second truth source is created.

## Layout

```
src/build_pair_dataset.py    # master-table builder (row = patient × interval × treatment)
src/normalization.py         # phase / treatment / endpoint / harmonization mappings
releases/v0.1.2/             # generated release (committed for review)
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
  referenced, not embedded.
- NeoTRIP is present only as an access-status record (see LABEL_RECOVERY_LEDGER.csv).

See `DATASET_CARD.md` for scope/limitations and `SCHEMA.md` for column definitions.
