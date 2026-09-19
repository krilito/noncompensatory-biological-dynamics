# Expression layer v0.3.1 — Level C1: technology-native quantitative expression

Level C1 answers one question per source: **what do these numbers mean, and what may
legally be done to them?** It then puts each source on canonical-gene coordinates.

It does **not** make sources comparable to each other. There is no cross-cohort
correction here: no ComBat, no global z-score, no joint quantile normalization, no
ranks, no paired deltas, no models. A value in one source is comparable only to other
values in that same source.

```text
source native matrix
  -> establish the scale from documentary provenance only   (SOURCE_EXPRESSION_CONTRACT.csv)
  -> validate that the declaration is not contradicted by the data
  -> apply the one transform that scale permits             (expression_transforms.py)
  -> feature -> canonical gene using the v0.2.1 feature maps
  -> collapse several features of one gene by the scale-appropriate rule
  -> source-specific canonical quantitative matrix          (local parquet, not committed)
```

### What v0.3.1 repaired

v0.3 computed counts as `mapped features -> sum to gene -> CPM -> log2(CPM + 1)`, where
the CPM denominator was the sum of the mapped features. That renormalized the retained
genes to exactly 1e6 per column and discarded the reads of every native feature that did
not enter the canonical gene space. v0.3.1 keeps the same route but takes the denominator
from the **complete native matrix**; the published invariant changed accordingly, and the
two checks are now reported separately ([Verification](#verification-that-nothing-was-normalized-twice)).

```text
source native matrix
  ├─ all native features                         -> native_count_mass   (the denominator)
  └─ mapped features -> sum duplicate per gene   -> / native_count_mass -> CPM -> log2(CPM+1)
```

## The rule book

`src/expression_transforms.py` is the mathematical authority. The routes are fixed:

| Declared scale | What is done | Gene collapse |
|---|---|---|
| `RAW_COUNTS` | sum duplicate features per gene, **then** CPM **on the native library denominator**, **then** `log2(CPM + 1)` | sum (counts are additive) |
| `TPM` | `log2(TPM + 1)` | median (abundance per feature is not additive across features) |
| `FPKM` | `log2(FPKM + 1)` | median |
| `LOG2_CPM`, `LOG2_TPM`, `LOG2_FPKM`, `LOG_EXPRESSION` | **identity** — never logged a second time | median |
| `ARRAY_NORMALIZED_LOG` | **identity** — already log2 (RMA, lumi) | median |
| `ARRAY_NORMALIZED_LINEAR` | `log2(x + 1)` — extension, ratified, see below | median |
| `LINEAR_ABUNDANCE_COMBAT_ADJUSTED` | `log2(x + 1)` — extension, ratified, see below | median |
| `UNKNOWN` | **nothing**; the source does not enter the layer | — |

Where several native features resolve to one canonical gene the collapse is the **median
of duplicate mapped native features per canonical gene**; for a gene-level source such as
GSE91061 (Entrez gene rows) that operation is the identity, since no gene has more than
one mapped feature. Calling it a "probe median" there would describe an event that does
not occur.

Two guards run before any transform: a declared count matrix must be non-negative and
at least 99% integer-like; a declared linear abundance or linear intensity must be
non-negative. A contradiction is reported as `SCALE_CONTRACT_FAILED`. The declared
scale is never silently rewritten to match the data, and numeric magnitude is never
used to *choose* a scale.

## Two declared-scale extensions, ratified by the Owner

The supplied vocabulary had no member for two real cases among the 23 sources. Both
additions are additive, both reuse one already-supplied transform, and neither performs
any renormalization. The Owner ratified both for v0.3.1:

* **`ARRAY_NORMALIZED_LINEAR`** — GSE20181 (`MAS5`) and GSE87455 (`quantile
  normalisation with BASE`) are documented as normalized array signals that are still
  *linear* (GSE87455 columns all sum to one fixed total). Sending them through
  `ARRAY_NORMALIZED_LOG` would have left raw intensities of up to 108,791 in a "log"
  matrix; converting them by library-size CPM would be a second normalization of data
  the platform already normalized. `log2(x + 1)` + median is the only route that is
  both legal and non-duplicating. **249 of the 688 ready pairs depend on it.**
* **`LINEAR_ABUNDANCE_COMBAT_ADJUSTED`** — GSE319641 (NeoTRIP, **150 pairs**) documents
  its own chain verbatim in GEO: *"TPM values were log2-transformed after adding an
  offset of 1, ComBat-corrected, and converted back to linear TPM."* The delivered file
  is therefore `2**y - 1` with `y = ComBat(log2(TPM + 1))`: it is no longer TPM (column
  sums are not 1e6, minimum is `-0.98`), so declaring `TPM` fails the non-negativity
  contract, and `log2(x + 1)` recovers exactly the author-computed `y`. That is not a
  normalization performed by us; it is the recovery of the log-space the author analysed.
  Validation rejects any value below the `-1` floor of that back-transform, and the
  source is flagged `author_batch_corrected=true`.

### Downstream red flag for the author-corrected source

`author_batch_corrected=true` is a modelling instruction, not a comment. Every per-source
metrics row now also carries two provenance fields:

| `preprocessing_scope` | `fold_isolation` | Meaning |
|---|---|---|
| `OUR_C1_PER_SAMPLE` | `FOLD_INDEPENDENT` | 22 sources: the value is a function of one sample column alone; nothing was fitted, so a train/test split cannot leak through it |
| `AUTHOR_FULL_SOURCE` | `NOT_ESTABLISHED` | `GSE319641_bulk`: the authors corrected batches **across the whole source**, so it may stay in the dataset as a real source representation but must never be presented as preprocessing we fitted inside a training fold |

If pre-ComBat TPM is located later, the primary model analysis should prefer it. Nothing
in C1 re-derives or undoes the author's correction.

## The RAW_COUNTS denominator: native, not mapped

A native feature that fails to map into the HGNC canonical space still consumed
sequencing depth. Two different quantities are therefore reported per source, and they
are not interchangeable:

* **feature mapping fraction** — what share of feature *IDs* became canonical genes;
* **canonical count mass fraction** — what share of *reads* those genes carried.

| Source | Feature mapping fraction | Canonical count mass fraction (median, range) | Effect of the v0.3 bug |
|---|---|---|---|
| GSE139533_bulk | 0.9995 | 0.9997 (0.9995–0.9998) | none worth naming |
| GSE310856_bulk | 0.6648 | 0.9923 (0.9909–0.9939) | genes inflated ≈0.8% |
| GSE152469_pseudobulk | 0.6679 | 0.9919 (0.9915–0.9922) | genes inflated ≈0.8% |
| GSE111014_pseudobulk | 0.6669 | 0.9928 (0.9871–0.9959) | genes inflated ≈0.7% |
| MGH_GSE115821_bulk | 0.2793 | **0.8452 (0.7577–0.9628)** | genes inflated ≈18%, up to 32% |
| MGH_GSE168204_bulk | 0.2797 | **0.9339 (0.6869–0.9700)** | genes inflated ≈7%, up to 46% |

MGH is the case that matters: ~72% of its 75,253 native rows are `NONHSAG*` IDs outside
the HGNC namespace, and those rows carry roughly 15% of the reads on a typical library
and up to a third on its worst. v0.3 deleted them from the denominator as well as from
the matrix, so a gene's linear CPM was too high by `1 / fraction`, i.e. by
`log2(1 / fraction)` on the written `log2(CPM + 1)` scale for any gene comfortably above
1 CPM (+0.24 bits at the MGH_GSE115821 median fraction, +0.54 bits at the worst library).
Every one of these numbers is recomputable from
`manifests/RAW_COUNT_LIBRARY_DENOMINATORS.csv.gz`, which lists `native_count_mass`,
`mapped_count_mass` and their ratio for all 108 count sample columns.

The worse consequence is not the offset but its variability: the fraction is a
per-library quantity (MGH_GSE168204 spans 0.687–0.970), so the old denominator injected
up to half a bit of **spurious between-sample difference** into each source, which is
exactly the signal a longitudinal delta would have picked up. Per-sample values are in
`manifests/RAW_COUNT_LIBRARY_DENOMINATORS.csv.gz`; no expected fraction is hard-coded
anywhere in the builder.


## Per-source adjudication

`Native scale` is the declared semantics; `what we did` is the transform actually
applied. `samples` are corpus sample rows bound to a matrix column.

| Source | Native scale (provenance) | What we did | Samples | Genes | Ready pairs |
|---|---|---|---|---|---|
| GSE91061_bulk | FPKM (file name) | `log2(FPKM+1)`, median | 107/108 | 21,813 | 43 |
| MORRISON_bulk | LOG2_CPM (file name: `...-logcpm-...`, "no batch correction") | **identity**, median | 442/442 | 21,050 | 71 |
| GSE139533_bulk | RAW_COUNTS (file name + `featureCounts`) | sum→native CPM→`log2(+1)` | 25/56 | 18,425 | 0 |
| GSE310856_bulk | RAW_COUNTS (file name + STAR/featureCounts) | sum→native CPM→`log2(+1)` | 33/33 | 39,362 | 22 |
| GSE319794_bulk | TPM (file name; every column sums to exactly 1e6) | `log2(TPM+1)`, median | 41/41 | 18,800 | 27 |
| GSE207422_bulk | LOG2_TPM (file name `...log2TPM`) | **identity**, median | 24/39 | 34,848 | 0 |
| GSE319641_bulk | COMBAT_ADJUSTED (GEO text, quoted above) | `log2(x+1)`, median | 401/401 | 23,884 | 150 |
| MGH_GSE115821_bulk | RAW_COUNTS (our HTSeq assembly) | sum→native CPM→`log2(+1)` | 13/21 | 20,997 | 6 |
| MGH_GSE168204_bulk | RAW_COUNTS (our HTSeq assembly) | sum→native CPM→`log2(+1)` | 22/24 | 21,023 | 9 |
| GIDE_bulk | **UNKNOWN** | **not processed** | 91/91 | — | 0 |
| GSE18728_array | ARRAY_NORMALIZED_LOG (RMA) | **identity**, median | 61/61 | 19,931 | 33 |
| GSE55374_array | ARRAY_NORMALIZED_LOG (lumi quantile) | **identity**, median | 36/36 | 20,550 | 23 |
| GSE20181_array | ARRAY_NORMALIZED_LINEAR (MAS5) | `log2(x+1)`, median | 176/176 | 12,515 | 114 |
| GSE87455_array | ARRAY_NORMALIZED_LINEAR (BASE quantile) | `log2(x+1)`, median | 275/275 | 20,550 | 135 |
| GSE3578_array | **UNKNOWN** | **not processed** | 0/78 | — | 0 |
| GSE65303_array | **UNKNOWN** | **not processed** | 18/18 | — | 0 |
| TRIO_US_B07_array | **UNKNOWN** (two-channel log ratios) | **not processed** | 1/1 | — | 0 |
| GSE111014_pseudobulk | RAW_COUNTS (our pseudobulk sums) | sum→native CPM→`log2(+1)` | 12/12 | 22,390 | 8 |
| GSE152469_pseudobulk | RAW_COUNTS (our UMI sums) | sum→native CPM→`log2(+1)` | 3/3 | 22,423 | 2 |
| GSE116256_pseudobulk | LOG2_CPM (`single_cell_qc.py` writes `log2(cpm+1)`) | **identity**, median | 35/36 | 25,850 | 19 |
| GSE123813_BCC_pseudobulk | LOG2_CPM (same) | **identity**, median | 22/22 | 18,763 | 11 |
| GSE123813_SCC_pseudobulk | LOG2_CPM (same) | **identity**, median | 8/8 | 15,629 | 4 |
| GSE165897_pseudobulk | LOG2_CPM (same) | **identity**, median | 22/22 | 23,154 | 11 |

Provenance for every row, quoted verbatim where it comes from GEO or a file name, is in
`SOURCE_EXPRESSION_CONTRACT.csv`.

### Why the four UNKNOWN sources are not processed

* **GIDE_bulk (16 pairs).** The file is titled `cancercell_normalized_counts_genenames`
  and is non-integer, while the only local documentation describes the accession as
  *"the read counts obtained from RNAseq data"* and names neither the normalization nor
  its divisor. `RAW_COUNTS` would fail the integer contract; any linear-abundance
  declaration would license a second, undocumented normalization of data whose name
  already says it was normalized. This is the source the double-normalization check
  exists for.
* **GSE3578_array.** Only *"global median normalization"* is documented for a 2004
  spotted cDNA/oligo platform; the delivered values (median 1.008, 0.47% negative,
  maximum 1221.5) are neither a documented linear intensity nor a documented log scale.
  It is also outside the canonical gene space (43.9% of features map).
* **GSE65303_array.** GEO names the software but never the scale of the delivered
  values. Inferring `log2` from a 1.4–23.1 range is exactly the forbidden shortcut.
* **TRIO_US_B07_array.** Agilent Feature Extraction with Linear/LOWESS plus Rosetta
  error-weighted averaging, centred on zero (50.5% negative, −2.85…2.40): a two-channel
  log-ratio quantity, so there is no single-sample abundance to place on a scale.

## Sample-to-column binding

C1 is where a MASTER sample id is first bound to a matrix column (v0.2 deliberately
declined to assert this). Each source declares one rule, and a binding is accepted only
when it is unique in both directions:

| Rule | Used by | Basis |
|---|---|---|
| `GSM_ID_HEADER` | all series-matrix sources | GEO accession is the column header |
| `NATIVE_ID_HEADER` | GSE91061, MORRISON, GSE319794, GSE207422, GSE139533, four pseudobulk sources | corpus native sample id is the column header |
| `GEO_SAMPLE_LABEL` | GSE310856, GSE319641 | the sample's own `!Sample_title` / `!Sample_description` value (`Library name: X` → `X`) |
| `MGH_LIBRARY_LABEL` | both MGH matrices | native label after dropping a trailing `.bam`/`_bam` and unifying `-`/`_` |
| `DERIVED_DONOR_DAY` | GSE111014_pseudobulk | our own generator's key `donor_D<days>` |

1,868 of 2,004 sample rows bind. The 136 that do not are all reported, never guessed:
120 are corpus rows naming several matrix columns
(`AMBIGUOUS_MULTIPLE_MATRIX_COLUMNS`) and 16 have no column at all in that file
(`NO_UNIQUE_MATRIX_COLUMN`, e.g. 15 GSE207422 immune-fraction samples that are absent
from the bulk matrix). See `manifests/SAMPLE_COLUMN_BINDING.csv.gz`. This costs **11
longitudinal pairs** (7 MGH_GSE115821, 2 MGH_GSE168204, 2 GSE139533).

### A corpus row naming several columns is not a set of technical replicates

Several libraries under one sample row can mean any of: one biospecimen sequenced
several times, several regions of one tumor, several pieces of one lesion, or several
assays of one visit. Collapsing the wrong one destroys the structure being measured —
GSE139533 is a multisampling study of intratumoral heterogeneity (128 tissue samples
from 44 tumors), so averaging its columns would erase the phenomenon it exists to
observe. C1 therefore collapses **nothing** and publishes
`manifests/SAMPLE_MULTI_COLUMN_ADJUDICATION.csv.gz`, one row per case, typed from
documentary evidence:

| `replicate_type` | Cases | What C1 may do |
|---|---|---|
| `TECHNICAL_LIBRARY_REPLICATE` | **0** | the only type eligible for a sum, and only for `RAW_COUNTS`, and only *before* CPM; the builder refuses to run if this type is ever assigned, because no such sum is implemented |
| `BIOLOGICAL_REGION_REPLICATE` | 31 (GSE139533_bulk) | preserve separately; never collapse |
| `BIOLOGICAL_ALIQUOT_REPLICATE` | 0 | preserve separately; never collapse |
| `DISTINCT_BIOSPECIMENS` | 0 | preserve separately; never collapse |
| `UNRESOLVED` | 89 (78 GSE3578_array, 8 MGH_GSE115821, 2 MGH_GSE168204, 1 GSE91061) | remain unbound |

Each row carries the patient uid where the interval table has one, the native timepoint,
the candidate columns, the GEO labels, the declared input scale, the corpus
`replicate_class`, and the quoted evidence (GEO `!Series_summary` / `!Sample_title` /
platform membership for GSE139533 and the MGH matrices; the absence of any aliquot
statement for GSE3578 and GSE91061). Losing 11 pairs to this rule is the intended
outcome: an unbound pair is recoverable by refining the corpus sample rows, a collapsed
one is not. For an already-normalized matrix (TPM/FPKM/logCPM/array) no generic mean or
sum rule is invented at all — such a case needs a source-provided aggregate or its
original raw counts.

## Verification that nothing was normalized twice

Two checks are published separately, because passing the first one says nothing about
the second. This is the distinction v0.3 blurred: all six counts sources passed with a
difference of exactly 0.0 while the denominator they were checked against was itself the
bug.

**`ROUTE_IMPLEMENTATION_CHECK`** — for each ready source the declared route is
recomputed *from the native matrix* and compared with the written matrix, restricted to
canonical genes served by exactly one mapped native feature, so neither the median nor
the sum collapse can hide a wrong transform. All 19 ready sources agree to within
**0.0** (columns `route_check_expected`, `route_check_genes`,
`route_check_max_abs_diff`); a disagreement sets `ROUTE_CHECK_FAILED` instead of
publishing. For counts sources the expected route string is now
`LOG2_OF_NATIVE_CPM_PLUS_1`.

**`LIBRARY_DENOMINATOR_CHECK`** — the deleted invariant was
`sum(2**x - 1) == 1e6` per column, which held by construction. What must hold instead,
for every RAW_COUNTS source and every sample column:

```text
sum(2**output - 1) == 1e6 * mapped_count_mass / native_count_mass
```

All six counts sources satisfy it, with a maximum relative deviation over all their
columns of `1.2e-13` (`canonical_cpm_mass_identity_holds`,
`library_denominator_max_rel_dev`); a violation sets `LIBRARY_DENOMINATOR_FAILED`. Their
column sums now read 757,709–999,807 instead of a suspiciously perfect 1,000,000.

The scale-independent part of the same evidence stays published for **all** ready
sources: `2**x - 1` summed per sample must not land on 1e6 for a source that arrived
already normalized, because that is what a second library-size normalization by us would
look like. The published `linear_library_sum_min/max` per ready source:

| Source class | Range of the per-column sum of `2**x - 1` |
|---|---|
| MORRISON log2 CPM | 0.61M–4.06M |
| GSE91061 FPKM | 0.31M–0.62M |
| GSE319641 author-ComBat linear TPM | 0.20M–4.76M |
| four author-log pseudobulk files | 0.89M–0.997M |
| microarray sources | 0.81M–8.57M |
| GSE319794 TPM | 0.993M–0.999M |
| RAW_COUNTS sources (ours, native denominator) | 0.687M–0.9998M |

Two of these sit close to a million, and for one of them that is expected: GSE319794 is
author-computed TPM, so its columns are near-per-million **because the author per-million
normalized them**, and we only logged them. The guard is an exact ±1 test on the sum,
which none of the 13 non-counts ready sources hits; had we renormalized any of them, they
would read as exactly 1e6.


## Headline

| Metric | Value |
|---|---|
| Sources declared | 23 |
| Sources quantitative-ready | **19** |
| Sources semantics-not-established (fail closed) | 4 |
| Sources with a contradicted declaration | 0 |
| Sources failing either published check | 0 |
| Samples in a ready matrix | 1,758 (of 2,004 declared) |
| Longitudinal pairs with t0 and t1 expression | 754 |
| — in the canonical gene space | 715 |
| **Pairs quantitative-ready** | **688** (492 patients) |
| — bulk RNA-seq / microarray / scRNA pseudobulk | 328 / 305 / 55 |
| — of those, clinical endpoint available | 394 |
| — of those, strict PR/CR-vs-PD eligible | 77 |
| — of those, frozen A/U/O eligible | 27 (all of them) |
| Multi-column sample cases adjudicated | 120, of which 0 collapsed |
| Lowest RAW_COUNTS canonical count mass fraction (source median) | 0.845 (MGH_GSE115821_bulk) |
| Derived matrix files | 19, 245.6 MiB, local only |

Rebuilding the six counts matrices changed no pair, sample or patient count: the fix
moves values, not membership.

Readiness is about expression, not labels. 150 of the 688 ready pairs are NeoTRIP
(GSE319641), whose clinical endpoint is `LABEL_NOT_FOUND` pending the access review
recorded in `LABEL_RECOVERY_LEDGER.csv`: they are quantitatively usable for
unsupervised work only. 443 of 2,311 native matrix columns belong to no corpus sample
and are dropped (they are recorded in the manifests, not silently discarded).

## Files

| File | Content |
|---|---|
| `SOURCE_EXPRESSION_CONTRACT.csv` | one row per source: modality, declared scale, provenance kind, verbatim provenance, probe-based, aggregation, binding rule, status |
| `SOURCE_EXPRESSION_METRICS.csv` | per source: features, columns, samples bound and multi-column cases, mapped/unmapped features, canonical genes, transform, input and output min/max, negative, zero and integer-like fractions, `linear_library_sum_min/max`, the native and canonical count mass fractions, `library_denominator_max_rel_dev`, `canonical_cpm_mass_identity_holds`, both checks, `preprocessing_scope`, `fold_isolation`, `author_batch_corrected`, output path and size, failure text |
| `EXPRESSION_LAYER_REPORT.json` | headline counts, per-class pair counts, scale vocabulary and its ratification, the two separated checks, canonical count-mass coverage per counts source, multi-column adjudication summary, downstream modelling provenance, double-log guards, redistribution policy, non-goals |
| `PAIR_EXPRESSION_COVERAGE.csv.gz` | one row per expression-bearing longitudinal pair: both endpoints, bound columns, quantitative status and reason |
| `manifests/SAMPLE_COLUMN_BINDING.csv.gz` | one row per sample: binding rule, native tokens, matched label, chosen column, binding status |
| `manifests/RAW_COUNT_LIBRARY_DENOMINATORS.csv.gz` | one row per RAW_COUNTS sample column (108 rows): native count mass, mapped count mass, canonical count mass fraction |
| `manifests/SAMPLE_MULTI_COLUMN_ADJUDICATION.csv.gz` | one row per multi-column case (120 rows): patient uid, timepoint, candidate columns, labels, input scale, corpus replicate class, quoted evidence, `replicate_type`, `recommended_action` |

## Where the matrices are, and why they are not committed

Real matrices are written to `longitudinal-data/data/processed/expression_v0.3/` as
`<source>.canonical_log_expression.parquet` (genes × samples, one file per ready
source, 245.6 MiB in total). That directory name was kept across the v0.3.1 repair and
its files were rewritten in place, so the matrices on disk are v0.3.1 matrices. That path
is ignored by git, and this repository carries no
per-accession licence ledger: GEO supplementary files are distributed under their
submitters' terms and a derived matrix inherits them. Redistribution therefore needs an
explicit Owner decision, so only code, contracts, metrics, manifests and documentation
are published.

## Rebuild

```bash
python src/build_gene_space.py        # v0.2.1 feature maps must exist first
python src/build_expression_layer.py  # rewrites this directory + local matrices
python -m pytest tests/test_expression_layer.py -q
```

## Non-goals

* No cross-source comparability: this layer is a precondition for Level C2
  (within-sample ranks), which is where a shared representation is first attempted.
* No normalization of anything already normalized, and no logarithm of anything already
  logged; both are checked per source in `SOURCE_EXPRESSION_METRICS.csv`.
* No renormalization of the canonical genes either: retaining a subset of a library's
  features does not make that subset the library, which is the v0.3 bug this release
  fixes.
* No collapsing of several matrix columns into one sample, for any reason, in C1.
* No choice of library-size, filter, gene universe, replicate collapse or model.
