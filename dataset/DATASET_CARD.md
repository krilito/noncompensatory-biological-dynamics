# PAIR Longitudinal Dataset Card — v0.1

## Summary

One row per patient × longitudinal sample interval × treatment context across 115
cohorts (17 cancer types) curated in the PAIR corpus. Generated 2026-09-19 from
`longitudinal-data/corpus` canonical tables + the patient-archive recovery pipeline.

## Contents

| Metric | Value |
|---|---|
| Master rows (longitudinal intervals) | 1111 |
| Resource patient entries | 2814 |
| Unique patients after confirmed dedup | 2684 |
| Paired patients (valid same-patient interval) | 909 |
| Endpoint-verified patients | 1122 |
| Paired + endpoint-verified patients | 553 |
| Strict A/U/O eligible patients (PR/CR-vs-PD, paired) | 58 unique patients / 99 resource entries |
| Frozen A/U/O paper subset (GSE91061) | 27 |
| Duplicate-resource groups | 131 |

## Provenance and identity

- `patient_uid` merges resource entries only on confirmed evidence: author subject IDs
  (GIDE/MORRISON_gide) or `global_patient_identity` confirmed same-patient roots
  (GSE91061/MORRISON_038). The legacy `GIDE:13` root, which wrongly merges the two
  author-separated subjects PD1_13 and ipiPD1_13, is deliberately not used.
- `duplicate_resource_group` distinguishes SAME_PATIENT groups (merged in `patient_uid`)
  from SAME_STUDY_RESOURCE groups (GSE20181/GSE5462, shared GSM accessions; patients
  kept separate pending a curation decision).

## Endpoints

Native semantics are preserved: CRS stays CRS (ordinal 1–3), pCR stays pCR, RECIST
stays RECIST, SD is a real label (STABLE in the additive `endpoint_harmonized` field),
GSE20181's ≥50% ultrasound-volume response keeps its own endpoint system. Recovered
labels (GSE165897 CRS, GSE179994 lesion-level RECIST, GSE120575 author RECIST binary,
GSE20181/GSE5462/GSE18728 semantics) are integrated; see LABEL_RECOVERY_LEDGER.csv.
Conflicting native values are kept as MULTIPLE_NATIVE_VALUES_REVIEW, never silently
resolved. No pseudo-labels were created.

## Limitations

- `delta_days` is NOT_MEASURED throughout: source corpora record order-only timing.
- Treatment attribution is interval-level (regimen before the t1 sample); several
  cohorts have unresolved per-patient regimens (NOT_MEASURED).
- Expression layers B–E (canonical gene space, within-cohort normalization, ranks,
  paired deltas) are not built in v0.1; the MASTER table references native
  expression files (Level A).
- NeoTRIP (251 patients) excluded pending access review.
- Strict A/U/O counts here are dataset-level; the frozen manuscript endpoint remains
  the GSE91061 27-patient subset.

## Ethics / access

Only redistributable public metadata and derived tables are included. No
controlled-access clinical material. NeoTRIP appears only as an access-status record.
