# PAIR Longitudinal Dataset Card — v0.1.2

## Summary

One row per patient × longitudinal sample interval × treatment context across 115
cohorts (17 cancer types) curated in the PAIR corpus. Generated 2026-09-19 from
`longitudinal-data/corpus` canonical tables + the patient-archive recovery pipeline.
v0.1.2 is a small identity repair on top of v0.1.1: GSE20181/GSE5462 patients
confirmed identical by shared GSM sample accessions are now merged at the
`patient_uid` level to prevent cross-resource train/test leakage; both
accession-specific resource records are retained. No new cohorts, no endpoint
changes.

## Contents

| Metric | Value |
|---|---|
| Master rows (longitudinal intervals) | 1111 |
| Resource patient entries | 2814 |
| Patient uids after confirmed same-patient merges | 2626 (NOT a final unique biological-patient count) |
| GSM-confirmed merged patient pairs (GSE20181/GSE5462) | 58 (resource records retained) |
| Paired patients (valid same-patient interval) | 909 |
| Endpoint-verified patients | 1122 |
| Paired + endpoint-verified patients | 553 |
| strict_prcr_vs_pd_eligible patients (cross-dataset) | 58 |
| frozen_auo_eligible patients (GSE91061 paper subset) | 27 |
| Interval endpoint status | bound: 738 rows; NOT_INTERVAL_RESOLVED: 13; LABEL_NOT_FOUND: 360 |

## v0.1.2 identity repair

The crosswalk (`gse20181_gse5462_crosswalk.csv`: 116 shared GSMs, 58 patient pairs,
native patient ids agree on all 116) confirms that the GSE20181 and GSE5462 resource
entries describe the same 58 biological patients. They now share one `patient_uid`
per pair (and one leakage group), so any train/test split over `patient_uid` cannot
place the two accessions of one patient on opposite sides. Both resource records
remain in the dataset (per-cohort rows, duplicate_resource_group =
`DUP_GSE20181_GSE5462`). All other v0.1.1 semantics and counts are unchanged
(patient uids 2684 → 2626).

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
   patient ids agree — merged at patient_uid level in v0.1.2).
5. **Treatment interval semantics.** `treatment_at_t1` is preserved; the derived
   `interval_treatment` exists only when a PRE_TREATMENT t0 supports the chronology
   (INTERVAL_ATTRIBUTED: 647 rows; otherwise INTERVAL_ATTRIBUTION_UNRESOLVED).

## Companion layer: v0.2 canonical gene space

`gene_space/` maps every feature of the 23 locally present expression files onto HGNC
approved genes, preserving each native identifier and recording ambiguity instead of
resolving it. Coverage is reported on two denominators (raw, and gene-addressable over
identifier classes present in the HGNC namespace); 22 of 23 sources pass the entry rule and
their shared core holds 9,338 genes over a 41,143-gene union. 715 of the 754 intervals that
have both endpoints locally available enter that gene space (496 patients). Metadata
semantics in this card are unchanged by it.

## Provenance and identity

- `patient_uid` merges resource entries only on confirmed evidence: author subject IDs
  (GIDE/MORRISON_gide), `global_patient_identity` confirmed same-patient roots
  (GSE91061/MORRISON_038), or shared GSM sample accessions with agreeing native
  patient ids (GSE20181/GSE5462). The legacy `GIDE:13` root, which wrongly merges the two
  author-separated subjects PD1_13 and ipiPD1_13, is deliberately not used.
- Recovered labels (GSE165897 CRS, GSE179994 lesion-level RECIST, GSE120575 author
  RECIST binary, GSE20181/GSE5462/GSE18728 semantics) are integrated; see
  LABEL_RECOVERY_LEDGER.csv. Conflicting native values are kept as
  MULTIPLE_NATIVE_VALUES_REVIEW. No pseudo-labels were created.

## Limitations

- `delta_days` is NOT_MEASURED throughout: source corpora record order-only timing.
- Interval endpoint binding resolves by t1-sample measurement or single-interval
  patients; 13 multi-interval rows remain NOT_INTERVAL_RESOLVED (fail-closed).
- Expression layers C–E (within-cohort normalization, ranks, paired deltas) are not
  built; the MASTER table references native expression files (Level A). Layer B, the
  canonical gene space, is delivered alongside v0.1.2 in `gene_space/` (see
  `gene_space/README.md`) and maps 715 of the 754 intervals that have both endpoints
  locally onto HGNC canonical genes. It identifies what each feature is; it does not make
  the measured values comparable.
- NeoTRIP (251 patients) excluded pending access review.
- `rebuild_snapshot/` contains the minimal redistributable inputs (public metadata
  only, 21 files, ~1.9 MB) needed to rebuild the release; see SNAPSHOT_MANIFEST.json.

## Ethics / access

Only redistributable public metadata and derived tables are included. No
controlled-access clinical material and no expression matrices. NeoTRIP appears
only as an access-status record.
