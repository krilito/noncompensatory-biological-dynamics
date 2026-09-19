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

## Companion layers

### v0.2.1 canonical gene space (level B)

`gene_space/` maps every feature of the 23 locally present expression files onto HGNC
approved genes, preserving each native identifier and recording ambiguity instead of
resolving it. Coverage is reported on two denominators (raw, and gene-addressable over
identifier classes present in the HGNC namespace); 22 of 23 sources pass the entry rule and
their shared core holds 9,338 genes over a 41,143-gene union. 715 of the 754 intervals that
have both endpoints locally available enter that gene space (496 patients). Metadata
semantics in this card are unchanged by it.

### v0.3.1 level C1: technology-native quantitative expression

`expression_layer/` gives each expression source a numerically valid scale *inside that
source* on canonical-gene coordinates, and is where a MASTER sample id is first bound to a
matrix column. Scales are adjudicated from documentary provenance only (file name, GEO
metadata, author scripts, our own derivation code) — never from numeric magnitude — and a
source whose semantics cannot be established fails closed. Counts are summed to gene, then
converted to CPM **on the denominator of the complete native matrix** — a feature that did
not map into the canonical gene space still consumed library depth, so the retained genes
are never renormalized to 1e6 (this is what v0.3.1 repaired); duplicate mapped features of
one gene take a median and are never summed; already-logarithmed matrices are never logged
again. 19 of 23 sources are quantitative-ready (1,758 samples), covering **688 of
the 754** longitudinal intervals (492 patients; 328 bulk RNA-seq / 305 microarray / 55
pseudobulk), of which 394 carry a clinical endpoint. A corpus sample row that names several
matrix columns is never collapsed: all 120 such cases are typed with quoted evidence in
`expression_layer/manifests/SAMPLE_MULTI_COLUMN_ADJUDICATION.csv.gz` (31 multi-sampling
regions, 89 unresolved), which is what costs 11 intervals. No cross-cohort correction is
performed: the matrices are explicitly **not** on one shared numerical scale, so values are
comparable only within a source. Level C2 is delivered alongside in `representation_layer/`
(see below); levels D and E remain unbuilt.
Derived matrices are local only (245.6 MiB) pending a redistribution decision.

### v0.3.2 pre-C2 engineering hardening (no change to C1 mathematics)

File locators in published artifacts are now **project-relative POSIX**, not machine
absolute paths: `samples.parquet::expression_file` and `endpoints.parquet::source_file`
lost 37,291 drive-qualified cells across the v0.1, v0.1.1 and v0.1.2 releases and the
committed `rebuild_snapshot/` corpus tables, with every other cell unchanged, and
`build_pair_dataset.portable_locator` keeps it that way on every rebuild. A locator that
cannot be expressed inside the project becomes `SOURCE_EXTERNAL` instead of an invented
path, so accession and table/column provenance stay the identity record. Public CI now also runs `dataset/tests`, whose guards were previously local-only,
and the GEO documenting-series lookup fails closed
(`MISSING_DOCUMENTING_SERIES_MATRIX` / `AMBIGUOUS_DOCUMENTING_SERIES_MATRIX`) instead of
silently returning no labels. Source counts, pair membership and every C1 value are
unchanged: 19 sources, 688 intervals, 492 patients.

### v0.4.1 level C2: two representations of the frozen C1 matrices

`representation_layer/` reads the 19 C1 matrices and never rebuilds them, and it ships two
things of different kinds. **C2-A is a fixed dataset representation**: every sample is ranked
against itself, `percentile = (average_rank - 0.5) / n_valid` (a midrank percentile, strictly
inside `(0, 1)`, ties averaged and never broken by gene order), over one shared coordinate
system — the frozen strict canonical core read from `gene_space/` (9,338 genes), not each
source's own gene list, because a `0.90` computed on 20,000 RNA-seq genes and one computed on
12,000 microarray genes are otherwise percentiles of different universes. No statistic is
shared between samples (measured per source: ranking one sample alone and inside its full
matrix agrees to 0.0), missing genes are never imputed and leave the denominator, and the
strict core turns out to be finite in essentially every ready sample (median coverage 1.0,
minimum 0.999036), so the coordinate really is shared. All 688 C1-ready intervals, 1,758
samples and 492 patients remain rank-ready: **C2 may only lose pairs, never recover one**, and
the builder raises if it ever does. **C2-B is a model-time contract, not a matrix**: training
median and 1.4826 × training MAD per gene, IQR/1.349 only when MAD is zero,
`UNUSABLE_CONSTANT_OR_SPARSE` when both are zero, and never an epsilon floor on the scale.
No whole-dataset standardized matrix exists anywhere in this repository, because fitting a
median or MAD over all samples fits it on the test patients. A held-out source with no
training sample may not have its scaler fitted from its own held-out samples either, and
GSE319641 — author ComBat over the whole cohort, `fold_isolation = NOT_ESTABLISHED` — is
refused by default unless a run declares itself a sensitivity/stress analysis; the flag
carries into all 150 of its rank-ready intervals. Rank resolution is source-specific and
published as such: a pseudobulk's zero block and an author's floored matrix are large ties
(786 distinct ranks in a median GSE116256 sample against 9,338 in GSE87455), inherited from
the source file rather than created here. Level C3 (paired state transitions) and levels D
and E remain unbuilt; derived rank matrices are local only (74.3 MiB).

**v0.4.1-C2 hardened the C2 API without touching its values**: the 19 rank matrices rebuild
byte-identically, which the newly published `rank_matrix_sha256` now proves. Three
second-order hazards found in code review are fail-closed. `upstream_fold_isolation` no
longer defaults to the safe answer, so a call that forgets to declare it raises instead of
quietly obtaining `FOLD_INDEPENDENT`. `apply_robust_standardizer` requires the caller to name
`expected_source_expression` and `expected_split_id`, and rejects cross-source, cross-split,
mixed and duplicated states, as well as a matrix with duplicate gene rows or duplicate
sample columns. The builder exits non-zero when any ready C1 source fails, after writing its
diagnostics. Two C2-A properties — every finite C1 value receives a rank, and every rank lies
strictly inside `(0, 1)` — are now build-stopping invariants rather than reported numbers. A
downstream level takes matrices from `SOURCE_RANK_METRICS.csv` via
`build_representation_layer.current_rank_matrices()` and never globs the matrix directory,
because a source that fails in a later build leaves an unmarked stale parquet behind.

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
- Expression levels C3 (paired deltas and state transitions), D and E are not built; the
  C2 rank representation is, and C2-B exists only as a fit/transform contract. Layer
  B, the canonical gene space, is delivered alongside v0.1.2 in `gene_space/` (see
  `gene_space/README.md`) and maps 715 of the 754 intervals that have both endpoints
  locally onto HGNC canonical genes: it identifies what each feature is. Layer C1, in
  `expression_layer/`, additionally gives each source a valid within-source numerical
  scale and binds samples to matrix columns, covering 688 of the 754 intervals; it does
  not make the measured values comparable across sources. Layer C2-A, in
  `representation_layer/`, puts all 688 of those intervals on one shared rank coordinate
  without making their absolute values comparable.
- NeoTRIP (251 patients) is excluded from clinical use pending access review: its
  expression is local and 150 of its intervals are quantitative-ready in C1, but every
  NeoTRIP interval endpoint is `LABEL_NOT_FOUND`.
- `rebuild_snapshot/` contains the minimal redistributable inputs (public metadata
  only, 21 files, ~1.9 MB) needed to rebuild the release; see SNAPSHOT_MANIFEST.json.

## Ethics / access

Only redistributable public metadata and derived tables are included. No
controlled-access clinical material and no expression matrices. NeoTRIP appears
only as an access-status record.
