# Representation layer v0.4.1-C2 — Level C2: two representations of one C1 matrix

Level C2 asks what may be done with the numbers C1 established. It ships **two
representations that are not the same kind of object**:

```text
                    C1 canonical quantitative expression  (v0.3.1, frozen)
                                     |
                    +----------------+----------------+
                    |                                 |
                 C2-A                               C2-B
        within-sample percentile rank        train-fold robust standardized state
        no statistic from any other sample   median + MAD of the TRAIN samples only
        a fixed dataset representation       a model-time fit/transform contract
        computed once, written as matrices   never computed once for all the data
                    |                                 |
        the shared cross-platform            the within-source quantitative geometry
             coordinate system                 that a fold of the model can use
```

The second line is the reason this layer is split at all. A robust standardization fitted
over every sample would have been fitted **using the test patients**, which is the leak the
whole release exists to prevent. So C2-B is delivered as code plus contract plus guards,
and no `standardized.parquet` exists anywhere in this repository.

`src/representation_transforms.py` is the mathematical authority. It was supplied by the
Owner and no rule in it may be swapped for `StandardScaler`, `scipy.stats.zscore`, a global
z-score, quantile normalization or ComBat. `src/build_representation_layer.py` applies the
C2-A half to the 19 `QUANTITATIVE_READY` C1 matrices, which it **reads and never rebuilds**.

### What v0.4.1 changed: where the guarantees live

Reviewing `2d22d59` found second-order API hazards that the mathematics did not cause and
that a future pipeline could trip over. The first pass added declarations to close them. The
second pass **deleted those declarations and made the same facts structural**, on the
principle that a constraint a caller can write is weaker than one the data flow enforces:

```text
C2-A   rank_one_sample(one sample's vector)  ->  no other sample is reachable from here
C2-B   RobustSplit(train ids, eval ids)      ->  the two sides cannot overlap, by construction
       RobustState(frame, source, split, fit_sample_ids, fit_gene_ids)
                                             ->  a scaler carries what it actually touched
       apply(matrix, state, split)           ->  leakage is shown by comparing memberships
```

| Hazard | How it is now prevented |
|---|---|
| a rank that depended on the cohort would have been caught by comparing one sample's ranks alone vs in-company, at a `1e-7` tolerance | the numeric tolerance and its column are **gone**. The matrix builder calls `rank_one_sample` once per column and that function's argument type is a single vector, so the leak is not measurable-away, it is unrepresentable. Tests delete, extremize and shuffle the other columns and assert the target column is bit-for-bit unchanged, and a spy asserts the primitive is never handed more than one sample |
| `upstream_fold_isolation` as a call argument, which a caller can type | `robust_split_for_source()` copies it from `../expression_layer/SOURCE_EXPRESSION_METRICS.csv`. A split that never loaded it carries an empty value, which fails closed at fit time, so nobody inherits the safe answer by forgetting or by guessing |
| a scaler's identity as two strings the caller re-declares at apply time | they live on the `RobustState` the fit returned, and `verify_fit_membership(state, split)` compares `fit_sample_ids` against `split.train_sample_ids` by id. Fitting is only permitted on a matrix whose columns **are** that membership, so a fit handed the full matrix — the bad implementation — raises on the spot |
| the builder recording a failed source, printing it, and returning 0 | artifacts are still written first, for diagnosis, and then the process exits non-zero; a build that did not complete publishes nothing, and `rank_matrices_for_build()` refuses it |

Membership everywhere is compared **by id, not by length**: deleting one gene and adding
another keeps a shape identical and is exactly the corruption a count check cannot see, so
`missing / unexpected / duplicated` id lists are what the guards compute and what their
error messages print.

**No value changed across either pass**: the 19 rank matrices rebuilt byte-for-byte identical
to `2d22d59`, and the published QC and coverage tables are content-identical. That equality
was verified while developing, by hashing the rebuilt files against a saved baseline — a
development check, deliberately not a runtime gate or a published field.

Two C2-A properties remain build-stopping rather than merely reported: a finite C1 gene value
that receives no rank, or a rank that touches 0 or 1, fails its source
(`rank_invariant_violation`), and a rank matrix whose gene rows or sample columns are not
exactly the frozen core and the C1-bound samples fails its source (`membership_check`). Both
run **before** a rank matrix is written, and both are tested against deliberately broken
rankers — a swapped gene row at unchanged shape, a dropped gene, a renamed sample — not just
against the good path.

### How a downstream level may find the rank matrices

`SOURCE_RANK_METRICS.csv` is the only truth about which files belong to this build. Use
`build_representation_layer.rank_matrices_for_build(build_id)`, which returns `build_id`,
`rank_matrix_relpath`, `samples_ranked` and the provenance flags per source; **never** glob
`expression_v0.4_rank/`. A source that fails in a later build leaves its previous parquet
sitting in that directory, and a glob would train on it. The loader therefore refuses a build
that is not the one the layer published, refuses one that did not complete, refuses one whose
source membership does not match C1's ready list, and never falls back to the files an
earlier build left behind. Nothing is deleted by a build: identity is the readable `build_id`
on each source row plus the `build_status` of the build in
`REPRESENTATION_LAYER_REPORT.json`, not a hash, so a consumer can state in prose which build
it trained on.

## C2-A — within-sample percentile over the frozen strict core

```text
percentile = (average_rank - 0.5) / n_valid         RANK_METHOD = average
```

`n_valid` is the number of finite genes **in that sample**, and only those genes are ranked.
The result is a midrank empirical percentile: strictly inside `(0, 1)`, with the extremes at
`0.5 / N` and `1 - 0.5 / N`, so nothing lands on 0 or 1 and no value is a hard ceiling.

The ranking universe is **one fixed gene list for every source**: the strict canonical core,
read at build time from `../gene_space/GENE_SPACE_GENES.csv.gz` where
`in_core_gene_space` is true — currently 9,338 genes. The number is never hard-coded.

That choice is the whole point of C2-A. If RNA-seq ranked 20,000 genes and a microarray
ranked 12,000, then a `0.90` in one source and a `0.90` in the other would be percentiles of
*different reference universes* and could not be placed on one axis. Over a shared 9,338-gene
coordinate, a rank means the same position on the same list everywhere. The core is used
**as a measurement coordinate, not as a modelling feature set** — feature selection is a
later decision that belongs inside a fold.

### Missing values are not imputed

A gene that the source does not measure stays `NaN` in the rank matrix, is excluded from the
sample's ranks, and shrinks that sample's `n_valid`. Nothing is filled with zero, and no
sample is dropped by a threshold invented here: coverage is reported per source, and the
builder's published `finite_fraction` is what a downstream analysis must read before
choosing anything.

Measured, and this is a good outcome: the strict core is essentially complete in every ready
source. 18 of 19 sources carry all 9,338 core genes as finite in all samples;
GSE111014_pseudobulk carries 9,329 (median `finite_fraction` 0.999036). So `min` over all
1,758 samples is 0.999036 and the median is 1.0 — the ranking universe is not merely
*nominated* to be shared, it is shared.

### Ties are averaged, and they matter here more than usual

Genes with equal values receive the same percentile, computed as their average rank; the
order in which genes happen to be listed can never break a tie. Two sources show why that
rule is not decorative:

| Source | median largest tie block | median distinct ranks per sample (of 9,338) |
|---|---|---|
| GSE87455_array, GSE319641_bulk | 0.0001 (one gene) | 9,338 — a full ordering |
| GSE55374_array | 0.0002 | 9,335 |
| GSE20181_array | 0.0012 | 6,315 |
| MORRISON_bulk | 0.0193 | 3,780 |
| GSE18728_array | 0.1418 | 7,877 |
| GSE116256_pseudobulk | 0.1132 | **786** |
| GSE123813_SCC_pseudobulk | 0.1213 | 2,403 |

The blocks are inherited, not created here: for the sources checked natively, the repeated
value is already in the author's own file — `0.466129` appears 1,773 times in the GSE116256
series file before any mapping (1,558 times after), and `-4.203801`, a distribution floor,
appears 4,507 times in MORRISON's (4,488 after). A pseudobulk's zero block and an author's
rounded or floored matrix are large tie blocks by construction. Breaking them by gene order
would manufacture an ordering that the measurement does not contain.

The consequence is honest and belongs in any downstream report: **C2-A's resolution is
source-specific.** A percentile vector that collapses 9,338 genes into 786 distinct values
carries far less ordering information than one that keeps all 9,338, which is why
`rank_unique_values` is published per sample and `largest_tie_block_fraction_median` per
source rather than being averaged away.

### What makes C2-A safe to freeze

Ranking one sample uses no information from another, and that is a property of the code
rather than a measured hope: `rank_one_sample` receives one `pd.Series`, and
`build_within_sample_rank_matrix` calls it once per column and assembles the results. There
is no source mean, source variance, cohort median or batch correction anywhere in the path,
and the tests delete the neighbouring columns, replace them with `±1e9`, and reorder them
without moving the target column by one bit. Ranks are also invariant to any strictly
monotone change of the underlying scale, which is what makes them indifferent to the C1 scale
differences between sources.

### Output

`<source>.canonical_core_rank.parquet`, genes × samples, `float32`, 19 files totalling
74.3 MiB, written to `longitudinal-data/data/processed/expression_v0.4_rank/`. That path is
ignored by git for the same reason as the C1 matrices: a derived matrix inherits the licence
terms of the GEO supplementary file it came from, and this repository carries no per-accession
ledger, so redistribution needs an explicit Owner decision.

## C2-B — a fit/transform contract, deliberately not a matrix

```text
split patients  ->  RobustSplit(source_expression, split_id,
                                train_sample_ids, eval_sample_ids,
                                upstream_fold_isolation  <- copied from C1 by
                                                           robust_split_for_source())
                    -> fit_robust_standardizer(matrix_of_those_train_samples, split)
                    -> RobustState(frame, source_expression, split_id,
                                  fit_sample_ids, fit_gene_ids)   # what it actually touched
                         +-> apply_robust_standardizer(matrix, state, split)
                             used for train, validation and test alike: the same state, the
                             same split, nothing estimated from what it transforms, and the
                             state's membership checked against the split's every time
```

* centre = **training median**; scale = **1.482602218505602 × training MAD**. Cancer
  transcriptomes have extreme values, so `mean/std` is not used.
* If MAD is zero but the gene still varies (`0 0 0 4` in the training fold is the shape of
  the problem), the fallback is **training IQR / 1.3489795003921634**, recorded as
  `scale_method = IQR_FALLBACK`.
* If both are zero the gene is `UNUSABLE_CONSTANT_OR_SPARSE` and stays `NaN` under
  transform. Scale is **never** floored at an epsilon: `max(scale, 1e-8)` would turn an
  almost constant gene into an enormous z-score, and the guard test asserts the resulting
  values stay ordinary.
* Too few finite training samples for a gene is `INSUFFICIENT_TRAIN_SAMPLES` (default
  minimum 4), never a silent fit.
* **A fit may only be handed its own training samples.** The matrix's columns must be
  exactly `split.train_sample_ids` — a missing id, a duplicated id or one extra column
  raises, so the full source matrix cannot be passed and the split treated as a comment.
* **A split cannot put one sample on both sides.** `RobustSplit` rejects an overlapping
  train/eval membership, a duplicated id list and an empty training side at construction,
  so `verify_fit_membership` comparing two lists is a complete leakage proof rather than a
  plausibility check.
* Genes a caller transforms must be exactly the genes that were fitted, and every sample
  column must belong to this split's train or eval side; duplicate gene rows and duplicate
  sample columns are rejected, since a repeated column weights one sample twice.

### The unseen-source rule

**C2-B is not a zero-shot cross-source representation.** If an entire source or cohort is
held out and contributes no training sample, its median and MAD must not be fitted from that
held-out source — that is the same leak in a smaller hat. Such an evaluation should use
C2-A ranks or another explicitly source-independent representation, and a paper must say so.

### GSE319641 (NeoTRIP) may not wear a fold-isolated label

Its C1 upstream is `preprocessing_scope = AUTHOR_FULL_SOURCE`, `fold_isolation =
NOT_ESTABLISHED`, `author_batch_corrected = true`: the authors ComBat-adjusted the whole
cohort, so the test patients were seen before this layer ever existed. Nobody types that
value at fit time — `robust_split_for_source('GSE319641_bulk', ...)` reads it out of the C1
manifest and puts it on the split, and `fit_robust_standardizer` then refuses the split with
a `ValueError`. A split built by hand that never loaded its isolation carries an empty value,
which is refused too: a missing declaration is never read as the safe one. Passing
`allow_nonisolated_upstream=True` is permitted only for an explicitly named sensitivity or
stress analysis, and its numbers may never be reported as a leakage-isolated main result.

The same warning travels with the C2-A ranks: 150 of the 688 rank-ready pairs come from
that source, and they carry `preprocessing_scope`, `fold_isolation` and
`author_batch_corrected` in every artifact of this layer. **Ranking is fold-independent;
this source's representation is not**, and a rank transformation cannot retroactively hide
what the author's ComBat already saw.

## Coverage

The starting set is the 688 C1 quantitative-ready intervals, and C2 may only *lose* pairs,
never recover one — the builder raises if a pair is rated ready that C1 did not.

| | C1 ready | C2-A rank ready |
|---|---|---|
| sources | 19 | 19 |
| samples | 1,758 | 1,758 |
| intervals (pairs) | 688 | **688** |
| patients | 492 | **492** |
| — bulk RNA-seq | 328 | 328 |
| — microarray | 305 | 305 |
| — scRNA-seq pseudobulk | 55 | 55 |

Nothing was lost, because nothing had to be: every ready C1 column is bound to exactly one
MASTER sample, every ready sample covers essentially the whole strict core, and a sample
with zero finite genes — the only way to fail — does not occur. `ENDPOINT_SAMPLE_NOT_RANKED`
is in the vocabulary and is empty in this build.

## Files

| File | Content |
|---|---|
| `README.md` | this document |
| `SOURCE_RANK_METRICS.csv` | one row per ready C1 source: universe size, samples ranked, core coverage min/median/max, tie-block and distinct-rank resolution, rank min/max and its `(0, 1)` flag, `float32` flag, the three provenance flags, `build_id`, local matrix path, size, failure text. This is the manifest a downstream level must load matrices from, and the only place a matrix locator is legitimate to read |
| `SAMPLE_RANK_QC.csv.gz` | one row per ranked sample: universe genes, finite genes, finite fraction, rank min/max, distinct rank values, provenance flags |
| `PAIR_RANK_COVERAGE.csv.gz` | one row per longitudinal interval: C1 status carried through, rank status and reason, provenance flags, endpoint eligibility |
| `REPRESENTATION_LAYER_REPORT.json` | headline counts, the C2-A definition, the C2-B contract including the unseen-source rule and the NeoTRIP guard, modality breakdown, per-source metrics, redistribution policy, non-goals |

## Rebuild

```bash
python src/build_expression_layer.py       # v0.3.1-C1 matrices must exist first
python src/build_representation_layer.py   # rewrites this directory + local rank matrices
python -m pytest tests/test_representation_layer.py -q
```

## Non-goals

* No paired delta, no Δrank, no t0/t1 state transition — that is level C3.
* No response model, no classification, no feature selection, no PCA.
* No batch correction, and no cross-cohort z-scoring, in either representation.
* No imputation of absent core genes, and no sample exclusion threshold.
* No C2-B matrix of any kind.
* Not cross-source comparability of measured values: C2-A buys a shared coordinate system,
  not a shared measurement scale, and a rank in one source is still not a rank-scaled value
  in another in any absolute sense.
