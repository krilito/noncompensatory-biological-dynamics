# Expression layer v0.3 — Level C1: technology-native quantitative expression

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

## The rule book

`src/expression_transforms.py` is the mathematical authority. The routes are fixed:

| Declared scale | What is done | Gene collapse |
|---|---|---|
| `RAW_COUNTS` | sum duplicate features per gene, **then** CPM, **then** `log2(CPM + 1)` | sum (counts are additive) |
| `TPM` | `log2(TPM + 1)` | median (abundance per feature is not additive across features) |
| `FPKM` | `log2(FPKM + 1)` | median |
| `LOG2_CPM`, `LOG2_TPM`, `LOG2_FPKM`, `LOG_EXPRESSION` | **identity** — never logged a second time | median |
| `ARRAY_NORMALIZED_LOG` | **identity** — already log2 (RMA, lumi) | median probe |
| `ARRAY_NORMALIZED_LINEAR` | `log2(x + 1)` — extension, see below | median probe |
| `LINEAR_ABUNDANCE_COMBAT_ADJUSTED` | `log2(x + 1)` — extension, see below | median |
| `UNKNOWN` | **nothing**; the source does not enter the layer | — |

Two guards run before any transform: a declared count matrix must be non-negative and
at least 99% integer-like; a declared linear abundance or linear intensity must be
non-negative. A contradiction is reported as `SCALE_CONTRACT_FAILED`. The declared
scale is never silently rewritten to match the data, and numeric magnitude is never
used to *choose* a scale.

## Two declared-scale extensions, flagged for architecture review

The supplied vocabulary had no member for two real cases among the 23 sources. Both
additions are additive, both reuse one already-supplied transform, and neither performs
any renormalization:

* **`ARRAY_NORMALIZED_LINEAR`** — GSE20181 (`MAS5`) and GSE87455 (`quantile
  normalisation with BASE`) are documented as normalized array signals that are still
  *linear* (GSE87455 columns all sum to one fixed total). Sending them through
  `ARRAY_NORMALIZED_LOG` would have left raw intensities of up to 108,791 in a "log"
  matrix; converting them by library-size CPM would be a second normalization of data
  the platform already normalized. `log2(x + 1)` + probe median is the only route that
  is both legal and non-duplicating. **249 of the 688 ready pairs depend on it.**
* **`LINEAR_ABUNDANCE_COMBAT_ADJUSTED`** — GSE319641 (NeoTRIP, **150 pairs**) documents
  its own chain verbatim in GEO: *"TPM values were log2-transformed after adding an
  offset of 1, ComBat-corrected, and converted back to linear TPM."* The delivered file
  is therefore `2**y - 1` with `y = ComBat(log2(TPM + 1))`: it is no longer TPM (column
  sums are not 1e6, minimum is `-0.98`), so declaring `TPM` fails the non-negativity
  contract, and `log2(x + 1)` recovers exactly the author-computed `y`. Validation
  rejects any value below the `-1` floor of that back-transform, and the source is
  flagged `author_batch_corrected=true` — **the batch correction is the author's, not
  ours, and this matrix must never be pooled as if it were uncorrected.**

## Per-source adjudication

`Native scale` is the declared semantics; `what we did` is the transform actually
applied. `samples` are corpus sample rows bound to a matrix column.

| Source | Native scale (provenance) | What we did | Samples | Genes | Ready pairs |
|---|---|---|---|---|---|
| GSE91061_bulk | FPKM (file name) | `log2(FPKM+1)`, median | 107/108 | 21,813 | 43 |
| MORRISON_bulk | LOG2_CPM (file name: `...-logcpm-...`, "no batch correction") | **identity**, median | 442/442 | 21,050 | 71 |
| GSE139533_bulk | RAW_COUNTS (file name + `featureCounts`) | sum→CPM→`log2(+1)` | 25/56 | 18,425 | 0 |
| GSE310856_bulk | RAW_COUNTS (file name + STAR/featureCounts) | sum→CPM→`log2(+1)` | 33/33 | 39,362 | 22 |
| GSE319794_bulk | TPM (file name; every column sums to exactly 1e6) | `log2(TPM+1)`, median | 41/41 | 18,800 | 27 |
| GSE207422_bulk | LOG2_TPM (file name `...log2TPM`) | **identity**, median | 24/39 | 34,848 | 0 |
| GSE319641_bulk | COMBAT_ADJUSTED (GEO text, quoted above) | `log2(x+1)`, median | 401/401 | 23,884 | 150 |
| MGH_GSE115821_bulk | RAW_COUNTS (our HTSeq assembly) | sum→CPM→`log2(+1)` | 13/21 | 20,997 | 6 |
| MGH_GSE168204_bulk | RAW_COUNTS (our HTSeq assembly) | sum→CPM→`log2(+1)` | 22/24 | 21,023 | 9 |
| GIDE_bulk | **UNKNOWN** | **not processed** | 91/91 | — | 0 |
| GSE18728_array | ARRAY_NORMALIZED_LOG (RMA) | **identity**, probe median | 61/61 | 19,931 | 33 |
| GSE55374_array | ARRAY_NORMALIZED_LOG (lumi quantile) | **identity**, probe median | 36/36 | 20,550 | 23 |
| GSE20181_array | ARRAY_NORMALIZED_LINEAR (MAS5) | `log2(x+1)`, probe median | 176/176 | 12,515 | 114 |
| GSE87455_array | ARRAY_NORMALIZED_LINEAR (BASE quantile) | `log2(x+1)`, probe median | 275/275 | 20,550 | 135 |
| GSE3578_array | **UNKNOWN** | **not processed** | 0/78 | — | 0 |
| GSE65303_array | **UNKNOWN** | **not processed** | 18/18 | — | 0 |
| TRIO_US_B07_array | **UNKNOWN** (two-channel log ratios) | **not processed** | 1/1 | — | 0 |
| GSE111014_pseudobulk | RAW_COUNTS (our pseudobulk sums) | sum→CPM→`log2(+1)` | 12/12 | 22,390 | 8 |
| GSE152469_pseudobulk | RAW_COUNTS (our UMI sums) | sum→CPM→`log2(+1)` | 3/3 | 22,423 | 2 |
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
120 are corpus rows that pool several sequencing libraries of one biological sample
(`AMBIGUOUS_MULTIPLE_MATRIX_COLUMNS`) — choosing or averaging among those libraries is a
modelling decision, so C1 refuses it — and 16 have no column at all in that file
(`NO_UNIQUE_MATRIX_COLUMN`, e.g. 15 GSE207422 immune-fraction samples that are absent
from the bulk matrix). See `manifests/SAMPLE_COLUMN_BINDING.csv.gz`. This costs **11
longitudinal pairs** (7 MGH_GSE115821, 2 MGH_GSE168204, 2 GSE139533).

## Verification that nothing was normalized twice

For each ready source the builder recomputes the declared route *from the native matrix*
and compares it with the written matrix, restricted to canonical genes served by exactly
one mapped native feature — so neither the median nor the sum collapse can hide a wrong
transform. All 19 ready sources agree to within 0.0 (columns `route_check_expected`,
`route_check_genes`, `route_check_max_abs_diff`), and a disagreement sets the source to
`ROUTE_CHECK_FAILED` rather than publishing it.

A second, scale-independent check: `2**x - 1` of the written matrix is summed per sample.
The six RAW_COUNTS sources come to exactly 1,000,000 per column (`cpm_invariant_holds`),
while every source that arrived already normalized comes to something else — MORRISON
0.61M–4.06M, the four author-log pseudobulk files 0.89M–0.996M, GSE91061 FPKM 0.31M–0.62M,
the arrays 2.7M–8.6M. Had any of them been re-normalized by us, they would have read as
exactly 1e6.

## Headline

| Metric | Value |
|---|---|
| Sources declared | 23 |
| Sources quantitative-ready | **19** |
| Sources semantics-not-established (fail closed) | 4 |
| Sources with a contradicted declaration | 0 |
| Samples in a ready matrix | 1,758 (of 2,004 declared) |
| Longitudinal pairs with t0 and t1 expression | 754 |
| — in the canonical gene space | 715 |
| **Pairs quantitative-ready** | **688** (492 patients) |
| — bulk RNA-seq / microarray / scRNA pseudobulk | 328 / 305 / 55 |
| — of those, clinical endpoint available | 394 |
| — of those, strict PR/CR-vs-PD eligible | 77 |
| — of those, frozen A/U/O eligible | 27 (all of them) |
| Derived matrix files | 19, 245.6 MiB, local only |

Readiness is about expression, not labels. 150 of the 688 ready pairs are NeoTRIP
(GSE319641), whose clinical endpoint is `LABEL_NOT_FOUND` pending the access review
recorded in `LABEL_RECOVERY_LEDGER.csv`: they are quantitatively usable for
unsupervised work only. 443 of 2,311 native matrix columns belong to no corpus sample
and are dropped (they are recorded in the manifests, not silently discarded).

## Files

| File | Content |
|---|---|
| `SOURCE_EXPRESSION_CONTRACT.csv` | one row per source: modality, declared scale, provenance kind, verbatim provenance, probe-based, aggregation, binding rule, status |
| `SOURCE_EXPRESSION_METRICS.csv` | per source: features, columns, samples bound, mapped/unmapped features, canonical genes, transform, input and output min/max, negative and zero fractions, integer-like fraction, output path and size, failure text |
| `EXPRESSION_LAYER_REPORT.json` | headline counts, per-class pair counts, scale vocabulary, the two extensions, double-log guards, redistribution policy, non-goals |
| `PAIR_EXPRESSION_COVERAGE.csv.gz` | one row per expression-bearing longitudinal pair: both endpoints, bound columns, quantitative status and reason |
| `manifests/SAMPLE_COLUMN_BINDING.csv.gz` | one row per sample: binding rule, native tokens, matched label, chosen column, binding status |

## Where the matrices are, and why they are not committed

Real matrices are written to `longitudinal-data/data/processed/expression_v0.3/` as
`<source>.canonical_log_expression.parquet` (genes × samples, one file per ready
source, 245.6 MiB in total). That path is ignored by git, and this repository carries no
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
* No choice of library-size, filter, gene universe, replicate collapse or model.
