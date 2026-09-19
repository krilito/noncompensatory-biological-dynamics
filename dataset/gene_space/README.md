# PAIR Canonical Gene Space — v0.2.1 (expression layer B)

This layer answers one question: **which biological feature does this expression column
represent?** It maps the heterogeneous feature identifiers of every locally available
expression file onto one coordinate system (HGNC approved gene symbols) while keeping the
native identifier of every feature.

It does **not** answer whether the measured numbers are comparable. No normalization, no
batch correction, no ranks and no paired deltas happen here — expression values stay in
their native units and technologies. Value comparability is the next layer.

## Inputs

| Input | Role | Provenance |
|---|---|---|
| `releases/v0.1.2/PAIR_LONGITUDINAL_MASTER.parquet` | longitudinal intervals + repository-relative expression references | committed |
| `releases/v0.1.2/samples.parquet` | binds each expression file to its cohort(s), modality and platform | committed |
| 23 locally present expression files | the features to map | referenced by path, never committed |
| `hgnc_complete_set.txt` | canonical gene table + symbol/alias/prev-symbol/Entrez/Ensembl crosswalks | `reference/hgnc_complete_set.txt` (16.9 MB, **not** committed; origin and SHA-256 recorded in `GENE_SPACE_REPORT.json`) |
| 6 GEO platform annotations (GPL96, GPL570, GPL2895, GPL6480, GPL10558, GPL16876) | probe → gene assignment | origin + SHA-256 recorded per platform |

## Artifacts

| File | Contents |
|---|---|
| `feature_map/<source>.feature_map.csv.gz` | one row per native feature of one expression source, with its canonical gene, its identifier class and the mapping route used |
| `SOURCE_MAPPING_METRICS.csv` | one row per expression source: both coverage denominators, status counts, unique canonical genes, entry-rule inputs and the identifier-class evidence |
| `GENE_SPACE_GENES.csv.gz` | one row per canonical gene in the union: HGNC/Ensembl/Entrez ids, locus group, how many sources / cohorts / modalities carry it, `in_core_gene_space` |
| `SOURCE_GENE_OVERLAP.csv` | source × source shared-gene counts — coverage evidence for choosing a modelling universe |
| `FEATURE_IDENTIFIER_INVENTORY.csv` | one row per cohort × modality slot in the corpus, with its expression files and their status |
| `PAIR_GENE_SPACE.csv.gz` | one row per longitudinal interval: its expression source and whether it enters the canonical gene space |
| `GENE_SPACE_REPORT.json` | headline metrics, reference digests, status vocabulary, low-coverage explanations, explicit non-goals |

## Feature-map columns

`cohort_code`, `source_expression`, `platform`, `original_feature_id`,
`original_feature_type`, `original_feature_annotation`, `feature_identifier_class`,
`gene_addressable`, `canonical_gene_id`, `hgnc_symbol`, `ensembl_gene_id`, `entrez_gene_id`,
`mapping_method`, `mapping_status`, `mapping_ambiguity`, `reference_version`.

`original_feature_id` and `original_feature_annotation` (the native gene symbols / gene ids
carried by the platform annotation) preserve the native representation; nothing is
overwritten, and a row that is called non-gene-addressable is still published with its
native identifier. `reference_version` names the HGNC snapshot and, for probe sources, the
GPL annotation snapshot by SHA-256 prefix.

## Status vocabulary

| Status | Meaning |
|---|---|
| `EXACT` | a direct identifier hit: approved HGNC symbol, stable Ensembl **gene** id, or Entrez gene id |
| `ALIAS` | resolved through an HGNC alias symbol |
| `PREVIOUS_SYMBOL` | resolved through an HGNC previous symbol |
| `PLATFORM_ANNOTATION` | a probe resolved through its GPL annotation to exactly one gene |
| `AMBIGUOUS` | the feature resolves to more than one canonical gene; `mapping_ambiguity` lists all candidates and **no** canonical gene is assigned |
| `UNMAPPED` | no canonical gene found; `mapping_method` records which route was tried and failed |

`mapping_method` carries the route (`ENSEMBL_GENE_ID`, `ENTREZ_GENE_ID`,
`AUTHOR_SYMBOL_FALLBACK:*`, `PLATFORM_ANNOTATION:GPL96:GENE_SYMBOL(previous_symbol)`,
`PLATFORM_ANNOTATION:GPL2895:ANNOTATION_WITHOUT_GENE_ASSIGNMENT`, …). Probes are never
mapped by expression correlation, genomic proximity or majority vote. Transcript accessions
(`ENST…`) are recorded as `TRANSCRIPT_ID_NOT_GENE_ID` and not mapped to a gene.

## Two coverage denominators, and which one decides entry

| Fraction | Denominator |
|---|---|
| `raw_mapping_fraction` | every native feature of the published matrix, as it stands |
| `gene_addressable_mapping_fraction` | only features whose native identifier belongs to a class represented in the HGNC gene namespace |

A class leaves the gene-addressable denominator only on **measured namespace evidence**: the
builder counts how many distinct tokens of that class appear anywhere in the HGNC namespace
(approved symbols + aliases + previous symbols) and requires the class to be present in
under 1% of its tokens, over at least 25 tokens. A class is never excluded because it
failed to map, so the denominator is not defined by the outcome it is used to judge.

Per-class evidence for every source is stored in `SOURCE_MAPPING_METRICS.csv`
(`non_addressable_classes`) and the class definitions in `GENE_SPACE_REPORT.json`.

**Entry rule.** A source enters the core when >=50% of its raw features map, **or** when
>=50% of its gene-addressable features map *and* the excluded foreign-namespace classes
account for >=90% of its unmapped rows. Both fractions are reported for every source either
way, and `raw_coverage_passes_entry_rule`,
`gene_addressable_coverage_passes_entry_rule` and `addressable_admission_evidenced` show
which route decided each source. Missing platform annotation is not an excuse: a probe
whose platform table carries no gene assignment stays inside the denominator, because the
missing information is the annotation and not a foreign namespace.

## Union and core, not weakest-platform intersection

* **Union** = every canonical gene found in at least one source: **41,143**.
* **Core** = genes found in *every source that passes the entry rule*: **9,338 genes**
  (9,330 protein-coding) over 22 of 23 sources.
* The naive intersection over all 23 sources is 8,841 — deliberately **not** used, so that
  one sparse 2004 array annotation cannot dictate the shared representation.
* Admitting the two MGH matrices on the gene-addressable denominator costs 16 genes of the
  core (9,354 → 9,338), which is the quantified answer to the concern that a weak source
  would otherwise drag the shared space down.
* Coverage tiers for the future modelling universe: 15,299 genes in >=90% of core-rule
  sources, 17,406 in >=75%, 20,218 in >=50%. **The 9,338-gene strict core is not the chosen
  model input**; `SOURCE_GENE_OVERLAP.csv` and the per-gene `n_cohorts` / `n_bulk_sources`
  / `n_array_sources` / `n_single_cell_sources` columns are the evidence for that later
  choice.

## Headline coverage

| Metric | Value |
|---|---|
| Cohorts with locally available expression | 27 (28 cohort × modality slots) |
| Expression sources (distinct files) | 23, all mapped |
| Sources passing the entry rule | 22 |
| Samples with a local expression file | 2,004 |
| Samples in the canonical gene space | 1,926 |
| Longitudinal intervals with t0 **and** t1 expression | 754 (535 patients) |
| Intervals usable through the canonical gene space | **715** (496 patients) |
| — bulk RNA-seq / microarray / scRNA pseudobulk | 355 / 305 / 55 |

## Source adjudications

* **MGH — feature universe RESOLVED, both matrices admitted.** The two files share one
  75,253-row gene universe in which 54,072 rows are named by the `NONHSAG<digits>` legacy
  non-coding assembly scheme. Of those tokens, **1** appears anywhere in the 103,218-token
  HGNC namespace, so the class is not a block of historical aliases we failed to resolve: it
  is a different naming namespace. Both matrices therefore keep a raw coverage of 27.9% /
  28.0% **and** report a gene-addressable coverage of 99.2% / 99.4% over 21,181 features,
  and enter the core on that basis. The 73 retired `LOC*` ids and 28 spreadsheet-dated
  symbols (`1-Mar`) were *not* excluded — they denote intended genes and count against the
  source. No row was deleted and no raw number was restated. MGH is still two sources, so
  the legacy annotation cannot hide the contemporary matrix.
* **GSE3578_array — stays outside the core (43.9% raw, 43.9% gene-addressable).** All
  54,359 probe ids are present in the GPL2895 table, but that 2004 spotted cDNA/oligo
  annotation assigns a gene identifier to only 25,552 of them (47.0%); the other 28,807
  are unannotated physical features, which the rule counts against the source. No mapping
  was invented for them: no expression correlation, no nearest gene, no majority vote.
* **GSE310856_bulk** shows the evidence gate working the other way: 14,337 of its features
  are clone-based temporary names present in the HGNC namespace in under 1% of tokens, so
  they leave the denominator (raw 66.5% → gene-addressable 87.5%). It already entered on
  raw coverage, so the admission did not depend on that exclusion.
* 33 intervals (snRNA-seq, spatial) have no local expression file at all and are recorded
  as `EXPRESSION_FILE_NOT_LOCAL`, not silently dropped.

## Rebuild

```bash
# the reference is not vendored: fetch it once and check the digest recorded in the report
curl -L -o gene_space/reference/hgnc_complete_set.txt \
  https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt
python src/build_gene_space.py     # rewrites gene_space/ from v0.1.2 + local expression files
python -m pytest tests/test_gene_space.py -q
```

The builder fails closed if a locally present expression file is not declared as a source,
or if a declared source has no local file, so a file can never be silently skipped.
