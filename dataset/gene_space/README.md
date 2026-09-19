# PAIR Canonical Gene Space — v0.2.0 (expression layer B)

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
| `feature_map/<source>.feature_map.csv.gz` | one row per native feature of one expression source, with its canonical gene and the mapping route used |
| `SOURCE_MAPPING_METRICS.csv` | one row per expression source: feature counts, status fractions, unique canonical genes, entry-rule result, unmapped-reason breakdown |
| `GENE_SPACE_GENES.csv.gz` | one row per canonical gene in the union: HGNC/Ensembl/Entrez ids, locus group, how many sources / cohorts / modalities carry it, `in_core_gene_space` |
| `SOURCE_GENE_OVERLAP.csv` | source × source shared-gene counts — coverage evidence for choosing a modelling universe |
| `FEATURE_IDENTIFIER_INVENTORY.csv` | one row per cohort × modality slot in the corpus, with its expression files and their status |
| `PAIR_GENE_SPACE.csv.gz` | one row per longitudinal interval: its expression source and whether it enters the canonical gene space |
| `GENE_SPACE_REPORT.json` | headline metrics, reference digests, status vocabulary, low-coverage explanations, explicit non-goals |

## Feature-map columns

`cohort_code`, `source_expression`, `platform`, `original_feature_id`,
`original_feature_type`, `original_feature_annotation`, `canonical_gene_id`, `hgnc_symbol`,
`ensembl_gene_id`, `entrez_gene_id`, `mapping_method`, `mapping_status`,
`mapping_ambiguity`, `reference_version`.

`original_feature_id` and `original_feature_annotation` (the native gene symbols / gene ids
carried by the platform annotation) preserve the native representation; nothing is
overwritten. `reference_version` names the HGNC snapshot and, for probe sources, the GPL
annotation snapshot by SHA-256 prefix.

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

## Union and core, not weakest-platform intersection

* **Union** = every canonical gene found in at least one source: **41,143**.
* **Core** = genes found in *every source that passes the entry rule* (>=50% of its
  features map to exactly one canonical gene): **9,354 genes** (9,345 protein-coding) over
  20 of 23 sources.
* The naive intersection over all 23 sources is 8,841 — deliberately **not** used, so that
  one sparse 2004 array annotation cannot dictate the shared representation.
* Coverage tiers for the future modelling universe: 15,322 genes in >=90% of core-rule
  sources, 17,589 in >=75%, 20,429 in >=50%. `SOURCE_GENE_OVERLAP.csv` and the per-gene
  `n_cohorts` / `n_bulk_sources` / `n_array_sources` / `n_single_cell_sources` columns are
  the evidence for choosing among them.

## Headline coverage

| Metric | Value |
|---|---|
| Cohorts with locally available expression | 27 (28 cohort × modality slots) |
| Expression sources (distinct files) | 23, all mapped |
| Sources passing the core entry rule | 20 |
| Samples with a local expression file | 2,004 |
| Samples in the canonical gene space | 1,881 |
| Longitudinal intervals with t0 **and** t1 expression | 754 (535 patients) |
| Intervals usable through the canonical gene space | **691** (481 patients) |
| — bulk RNA-seq / microarray / scRNA pseudobulk | 331 / 305 / 55 |

## Low-coverage sources are kept and explained

None of these is hidden: each keeps its own full feature map and contributes to the union;
they are excluded from the *core* only.

* **GSE3578_array** (GPL2895, 43.9% mapped): all 54,359 probe ids are present in the
  platform table, but that 2004 spotted cDNA/oligo annotation assigns a gene identifier to
  only 25,552 of them (47.0%). The ceiling is the annotation, not the crosswalk.
* **MGH_GSE115821_bulk** (27.9%) and **MGH_GSE168204_bulk** (28.0%): both files share one
  legacy HTSeq gene universe in which 54,071 rows are `NONHSAG*` non-coding identifiers
  that were never given an HGNC-approved symbol (plus 73 retired `LOC*` ids and 28
  Excel-date-corrupted symbols such as `1-Mar`). Excluding those tokens the same files map
  at 99.2% and 99.4%. MGH therefore appears as **two** sources, so a legacy annotation
  cannot hide the contemporary matrix.
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
