# PAIR Longitudinal Dataset Card — v0.1.1

## Summary

One row per patient × longitudinal sample interval × treatment context across 115
cohorts (17 cancer types) curated in the PAIR corpus. Generated 2026-09-19 from
`longitudinal-data/corpus` canonical tables + the patient-archive recovery pipeline.
v0.1.1 is a semantic repair of v0.1 (endpoint-to-interval binding, endpoint
semantics separation, eligibility and patient-count metric clarification); no new
cohorts were added.

## Contents

| Metric | Value |
|---|---|
| Master rows (longitudinal intervals) | 1111 |
| Resource patient entries | 2814 |
| Patient uids after confirmed same-patient merges | 2684 (NOT a final unique biological-patient count) |
| Unresolved same-study duplicate groups | 1 (GSE20181/GSE5462; 58 overlapping patients by shared GSM, records not merged) |
| Paired patients (valid same-patient interval) | 909 |
| Endpoint-verified patients | 1122 |
| Paired + endpoint-verified patients | 553 |
| strict_prcr_vs_pd_eligible patients (cross-dataset) | 58 |
| frozen_auo_eligible patients (GSE91061 paper subset) | 27 |
| Interval endpoint status | bound: 738 rows; NOT_INTERVAL_RESOLVED: 13; LABEL_NOT_FOUND: 360 |

## v0.1.1 semantic repairs

1. **Endpoint-to-interval binding.** Patient-level outcomes are kept in
   `patient_endpoint_*` fields. An interval receives `interval_endpoint_*` values
   only when the endpoint was measured on its t1 sample or the patient has exactly
   one valid interval; otherwise `interval_endpoint_status = NOT_INTERVAL_RESOLVED`
   and nothing is copied across transitions.
2. **Eligibility fields.** `strict_prcr_vs_pd_eligible` (cross-dataset, interval-bound
   PR/CR-vs-PD; 58 patients) is separate from `frozen_auo_eligible` (the frozen
   GSE91061 27-patient paper subset; the frozen endpoint itself is unchanged).
3. **Endpoint semantics.** Explicit `family_native / family_harmonized / system /
   native / harmonized` fields at both patient and interval level; native values
   are never replaced. GSE120575's author dichotomy is classified
   `RECIST_BINARY_AUTHOR_DICHOTOMY` (RECIST-derived, PMC6641984), not a generic
   author composite.
4. **Patient-count terminology.** `patient_uids_after_confirmed_same_patient_merges`
   replaces "unique patients"; unresolved same-study duplicate resources are
   reported separately with an explicit GSM-linked sample crosswalk
   (`gse20181_gse5462_crosswalk.csv`: 116 shared GSMs, 58 patient pairs, all native
   patient ids agree — merge remains a curation decision, not performed).
5. **Treatment interval semantics.** `treatment_at_t1` is preserved; the derived
   `interval_treatment` exists only when a PRE_TREATMENT t0 supports the chronology
   (INTERVAL_ATTRIBUTED: 647 rows; otherwise INTERVAL_ATTRIBUTION_UNRESOLVED).

## Provenance and identity

- `patient_uid` merges resource entries only on confirmed evidence: author subject IDs
  (GIDE/MORRISON_gide) or `global_patient_identity` confirmed same-patient roots
  (GSE91061/MORRISON_038). The legacy `GIDE:13` root, which wrongly merges the two
  author-separated subjects PD1_13 and ipiPD1_13, is deliberately not used.
- Recovered labels (GSE165897 CRS, GSE179994 lesion-level RECIST, GSE120575 author
  RECIST binary, GSE20181/GSE5462/GSE18728 semantics) are integrated; see
  LABEL_RECOVERY_LEDGER.csv. Conflicting native values are kept as
  MULTIPLE_NATIVE_VALUES_REVIEW. No pseudo-labels were created.

## Limitations

- `delta_days` is NOT_MEASURED throughout: source corpora record order-only timing.
- Interval endpoint binding resolves by t1-sample measurement or single-interval
  patients; 13 multi-interval rows remain NOT_INTERVAL_RESOLVED (fail-closed).
- Expression layers B–E (canonical gene space, within-cohort normalization, ranks,
  paired deltas) are not built in v0.1.x; the MASTER table references native
  expression files (Level A).
- NeoTRIP (251 patients) excluded pending access review.
- `rebuild_snapshot/` contains the minimal redistributable inputs (public metadata
  only, 21 files, ~1.9 MB) needed to rebuild the release; see SNAPSHOT_MANIFEST.json.

## Ethics / access

Only redistributable public metadata and derived tables are included. No
controlled-access clinical material and no expression matrices. NeoTRIP appears
only as an access-status record.
