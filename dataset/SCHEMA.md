# PAIR_LONGITUDINAL_MASTER Schema — v0.1.2

Row = one patient × one longitudinal sample interval × one treatment context.
`pair_uid` equals the corpus `transition_id` and is unique.

## Identity
| Column | Description |
|---|---|
| pair_uid | unique row id (= corpus transition_id) |
| cohort_code | corpus cohort code |
| patient_uid | stable internal id; merges resource entries only on confirmed identity evidence (author subject IDs, confirmed global_patient_identity roots, or shared GSM accessions with agreeing native patient ids — the GSE20181/GSE5462 merge added in v0.1.2) |
| native_patient_id | author-native patient id (author-corrected for GIDE/MORRISON_gide) |
| identity_status | CORPUS_NATIVE_ID / AUTHOR_SUBJECT_ID / SOURCE_PATIENT_ID_UNVERIFIED / ... |
| duplicate_resource_group | non-empty when the entry belongs to a confirmed duplicate-resource group |

## Disease
| Column | Description |
|---|---|
| cancer_type_native | native cancer label |
| cancer_type | harmonized cancer type |
| disease_subtype | OncoTree code where available |

## Longitudinal samples
| Column | Description |
|---|---|
| sample_t0 / sample_t1 | corpus sample ids (t0 before t1) |
| time_t0_native / time_t1_native | native timing strings, preserved verbatim |
| phase_t0 / phase_t1 | normalized phases: PRE_TREATMENT, ON_TREATMENT_EARLY, ON_TREATMENT_LATE, POST_TREATMENT, SURGERY, PROGRESSION, RECURRENCE, UNKNOWN |
| delta_days | exact interval in days; NOT_MEASURED everywhere in v0.1.x (sources are order-only) |

## Treatment
| Column | Description |
|---|---|
| treatment_at_t1 | native regimen recorded for the t1 sample ("treatment before t1") |
| interval_treatment | native regimen attributed to this interval; only when t0 is PRE_TREATMENT and t1 is on-treatment (treatment started after t0) |
| interval_treatment_status | INTERVAL_ATTRIBUTED / INTERVAL_ATTRIBUTION_UNRESOLVED |
| drug | parsed drug names (lexicon), else NOT_PARSED |
| drug_class | normalized classes, e.g. anti_PD1 + anti_CTLA4 |
| treatment_family | ICB, DUAL_ICB, ICB_CHEMOTHERAPY, ICB_COMBINATION, CHEMOTHERAPY, ENDOCRINE_THERAPY, TARGETED_THERAPY, RADIOTHERAPY, CHEMORADIATION, CELLULAR_THERAPY, OTHER, NOT_MEASURED |

## Clinical endpoint — patient level (never dropped)
Patient-level summary over all endpoint records linked to the patient. Kept even
when interval attribution is impossible.
| Column | Description |
|---|---|
| patient_endpoint_family_native | corpus endpoint families, pipe-joined (RECIST_RESPONSE, PCR, CHEMOTHERAPY_RESPONSE_SCORE, LESION_LEVEL_RECIST_RESPONSE, RECIST_BINARY_RESPONSE, CLINICAL_RESPONSE_ULTRASOUND_VOLUME, CLINICAL_RESPONSE_COMPOSITE_NEOADJUVANT, ENDPOINT_SEMANTICS_UNCERTAIN) |
| patient_endpoint_family_harmonized | harmonized families: radiologic_response, pathological_response, clinical_composite_response, recurrence, other_clinical_endpoint |
| patient_endpoint_system | RECIST, RECIST_1_1, RECIST_BINARY_AUTHOR_DICHOTOMY, pCR, CRS, ultrasound_volume_response, author_defined_composite, recurrence, uncertain |
| patient_endpoint_native | native label values, pipe-joined; conflicts visible, never collapsed |
| patient_endpoint_harmonized | additive only: CR/PR→RESPONSE, SD→STABLE, PD→PROGRESSION for standard RECIST_RESPONSE; NOT_HARMONIZED otherwise; never replaces the native label |
| patient_endpoint_assessment_time | assessment timing labels |
| patient_endpoint_status | fail-closed adjudication: OBSERVED_STRICT_RECIST, OBSERVED_PCR, OBSERVED_OTHER_ENDPOINT, MULTIPLE_NATIVE_VALUES_REVIEW, ENDPOINT_DEFINITION_UNRESOLVED, LABEL_NOT_FOUND, MAPPING_NOT_FOUND, PUBLIC_FILE_ACCESS_CONFLICT |
| patient_endpoint_verification_status | how endpoint rows joined to the patient (EXACT_SAMPLE:n / EXACT_PATIENT:n / NOT_LINKED) |

## Clinical endpoint — interval level (v0.1.1)
An interval receives `interval_endpoint_*` values ONLY when the endpoint can be
attributed to that interval: it was measured on the t1 sample of the interval, or
the patient has exactly one valid interval (unambiguous). Otherwise the outcome
stays patient-level and the interval is marked `NOT_INTERVAL_RESOLVED`; the
patient summary is never copied onto all transitions.
| Column | Description |
|---|---|
| interval_endpoint_status | same vocabulary as patient status, plus NOT_INTERVAL_RESOLVED; empty interval value fields whenever unresolved |
| interval_endpoint_family_native / _family_harmonized / _system | same semantics as the patient-level fields, restricted to bound endpoint rows |
| interval_endpoint_native | native values of the bound endpoint rows |
| interval_endpoint_harmonized | additive harmonization of the bound rows |
| interval_endpoint_assessment_time | timing labels of the bound rows |

## Expression
| Column | Description |
|---|---|
| expression_modality | bulk RNA-seq / microarray / two-channel microarray / scRNA-seq / snRNA-seq / spatial |
| platform | GEO platform id or descriptor |
| expression_t0_ref / expression_t1_ref | repository-relative path to the canonical expression file (Level A; no global normalization) |
| gene_space | PROBE_BASED / GENE_LEVEL / GENE_LEVEL_CELLULAR / GENE_LEVEL_SPATIAL |

## Provenance
| Column | Description |
|---|---|
| source_accession | dataset id (GSE… / cohort code) |
| source_publication | DOI/PMID for cohorts with a verified primary-paper definition, else empty |
| source_table | corpus source table of the row |
| provenance_id | joins to provenance.parquet (per-row source files, recovery flags) |

## Analysis views
| Column | Description |
|---|---|
| paired_valid | t0 and t1 samples map to the same verified patient |
| clinical_endpoint_available | patient-level endpoint_status starts with OBSERVED_ |
| strict_prcr_vs_pd_eligible | paired_valid AND the interval-bound endpoint status is OBSERVED_STRICT_RECIST; allowed across cohorts (58 unique patients in v0.1.2) |
| frozen_auo_eligible | the frozen GSE91061 paper subset only: strict interval-bound PR/CR-vs-PD in GSE91061 (27 unique patients, matching the frozen manuscript endpoint) |

No analysis-view column exists yet. When views are introduced they must not restamp
`cancer_type`: PAIR is a longitudinal cancer atlas, so a pan-cancer view includes cohorts
such as GSE139533 (glioblastoma, `cancer_type = Glioblastoma`) and a melanoma-only view
excludes them by filter, never by editing the record.

## Path locators and portability
| Column | Locator form |
|---|---|
| MASTER `expression_t0_ref` / `expression_t1_ref` | project-relative POSIX, or `NOT_AVAILABLE` |
| `samples.parquet::expression_file` | project-relative POSIX, or `NOT_FOUND` |
| `endpoints.parquet::source_file` | project-relative POSIX, or `SOURCE_EXTERNAL` |
| `gene_space/**`, `expression_layer/**` | project-relative POSIX only |

A published artifact says where a record's source lives *inside this project*; it never says
where one person's machine keeps it. `build_pair_dataset.portable_locator` performs the
conversion and `portable_frame` applies it to a redistributed snapshot; a value that is
already relative is read as project-relative and is never re-resolved against the working
directory, and an absolute path outside the project tree becomes `SOURCE_EXTERNAL` instead of
being rewritten into an invented path. Presence is a separate fact, recorded in
`source_path_status` (`LOCAL_SOURCE_PRESENT` / `SOURCE_PATH_UNAVAILABLE`); identity stays in
`cohort_code`, `source_accession`, `source_table`, `source_column` and `source_url`.
`dataset/tests/test_master_invariants.py` fails the build if any published tabular
artifact carries a drive-qualified, UNC or personal-directory path. In v0.3.2 the locator
columns of the already-published v0.1 and v0.1.1 releases and of `rebuild_snapshot/` were
rewritten in place; every other cell of those tables is unchanged.

## Release files
`patients.parquet` (resource-level patient entries), `samples.parquet`,
`endpoints.parquet` (linked endpoint rows incl. recovery sources),
`provenance.parquet`, `duplicate_groups.csv`,
`gse20181_gse5462_crosswalk.csv` (sample-level crosswalk for the same-study
duplicate resources, linked by shared GSM accessions only; v0.1.2 merges the 58
GSM-confirmed patient pairs at patient_uid level while retaining both resource
records),
`LABEL_RECOVERY_LEDGER.csv`, `DATASET_SUMMARY.json`,
`rebuild_snapshot/` (minimal redistributable input files + SNAPSHOT_MANIFEST.json
with public origins, sufficient to rebuild the release without the private workspace).

## Companion layer files

`gene_space/` (v0.2.1 canonical gene space, level B) is documented in
`gene_space/README.md`: its feature maps, source metrics, gene-level union/core table,
source overlap matrix, cohort inventory and per-interval `PAIR_GENE_SPACE.csv.gz`.

`expression_layer/` (v0.3.1 level C1) is documented in `expression_layer/README.md`:
`SOURCE_EXPRESSION_CONTRACT.csv` (declared scale plus verbatim provenance per source),
`SOURCE_EXPRESSION_METRICS.csv` (transform, input and output value ranges, feature and
sample counts, native and canonical count mass per source, both published checks,
downstream preprocessing provenance, matrix size, failure text),
`EXPRESSION_LAYER_REPORT.json`,
`PAIR_EXPRESSION_COVERAGE.csv.gz` (per-interval quantitative status),
`manifests/SAMPLE_COLUMN_BINDING.csv.gz` (one row per sample: binding rule, native
tokens, matched label, chosen matrix column, binding status),
`manifests/RAW_COUNT_LIBRARY_DENOMINATORS.csv.gz` (one row per count sample column: native
and mapped count mass and their ratio) and
`manifests/SAMPLE_MULTI_COLUMN_ADJUDICATION.csv.gz` (one row per sample that names several
matrix columns: evidence, replicate type, recommended action).

`representation_layer/` (v0.4.1 level C2) is documented in `representation_layer/README.md`:
`SOURCE_RANK_METRICS.csv` (per source: ranking universe, samples ranked, strict-core
coverage, tie-block and distinct-rank resolution, rank bounds and the `(0, 1)` flag, the three
preprocessing-provenance flags, and `build_id` plus the local matrix locator, size and failure
text — this table is the manifest a downstream level must load rank matrices from, through
`rank_matrices_for_build(build_id)`, never a glob of the matrix directory),
`SAMPLE_RANK_QC.csv.gz` (per sample: finite core genes, finite fraction, rank min/max,
distinct rank values), `PAIR_RANK_COVERAGE.csv.gz` (per interval: C1 status carried through,
rank status and reason, provenance flags) and `REPRESENTATION_LAYER_REPORT.json` (the C2-A
definition and the C2-B contract, including the unseen-source rule, the NeoTRIP
non-isolation guard and the matrix-consumption rule).

No companion layer adds columns to the MASTER schema or changes any v0.1.x semantic field.
Level C1 introduces the status vocabulary `QUANTITATIVE_READY`,
`SEMANTICS_NOT_ESTABLISHED`, `SCALE_CONTRACT_FAILED`, `SAMPLE_BINDING_FAILED`,
`ROUTE_CHECK_FAILED`, `LIBRARY_DENOMINATOR_FAILED`,
`SOURCE_NOT_RESOLVED_IN_GENE_SPACE` and `ENDPOINT_SAMPLE_NOT_BOUND_TO_MATRIX_COLUMN`,
and the declared-scale vocabulary in `expression_layer/README.md`. Level C2 introduces
`RANK_READY`, `ENDPOINT_SAMPLE_NOT_RANKED` and `NOT_C1_QUANTITATIVE_READY`, the C2-B
`fit_status` vocabulary `READY`, `INSUFFICIENT_TRAIN_SAMPLES`,
`UNUSABLE_CONSTANT_OR_SPARSE` and the `scale_method` vocabulary `MAD`, `IQR_FALLBACK`,
`UNUSABLE`, and adds no column to any earlier artifact. C2's build vocabulary is
`COMPLETE` / `INCOMPLETE` in `build_status`, published per source row and in the report
headline: a build in which every ready C1 source ranked is `COMPLETE`, and only a `COMPLETE`
build is served to a consumer.
