# PAIR_LONGITUDINAL_MASTER Schema — v0.1

Row = one patient × one longitudinal sample interval × one treatment context.
`pair_uid` equals the corpus `transition_id` and is unique.

## Identity
| Column | Description |
|---|---|
| pair_uid | unique row id (= corpus transition_id) |
| cohort_code | corpus cohort code |
| patient_uid | stable internal id; merges resource entries only on confirmed identity evidence |
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
| delta_days | exact interval in days; NOT_MEASURED everywhere in v0.1 (sources are order-only) |

## Treatment
| Column | Description |
|---|---|
| treatment_native | native regimen string for the interval (regimen before the t1 sample) |
| drug | parsed drug names (lexicon), else NOT_PARSED |
| drug_class | normalized classes, e.g. anti_PD1 + anti_CTLA4 |
| treatment_family | ICB, DUAL_ICB, ICB_CHEMOTHERAPY, ICB_COMBINATION, CHEMOTHERAPY, ENDOCRINE_THERAPY, TARGETED_THERAPY, RADIOTHERAPY, CHEMORADIATION, CELLULAR_THERAPY, OTHER, NOT_MEASURED |

## Clinical endpoint (patient-level summary projected onto each interval)
| Column | Description |
|---|---|
| endpoint_family | corpus endpoint families, pipe-joined (RECIST_RESPONSE, PCR, CHEMOTHERAPY_RESPONSE_SCORE, LESION_LEVEL_RECIST_RESPONSE, RECIST_BINARY_RESPONSE, CLINICAL_RESPONSE_ULTRASOUND_VOLUME, CLINICAL_RESPONSE_COMPOSITE_NEOADJUVANT, ENDPOINT_SEMANTICS_UNCERTAIN) |
| endpoint_system | RECIST, RECIST_1_1, pCR, CRS, ultrasound_volume_response, author_defined_composite, recurrence, uncertain |
| endpoint_native | native label values, pipe-joined; conflicts visible, never collapsed |
| endpoint_harmonized | additive only: CR/PR→RESPONSE, SD→STABLE, PD→PROGRESSION for standard RECIST_RESPONSE; NOT_HARMONIZED otherwise; never replaces the native label |
| endpoint_assessment_time | assessment timing labels from the endpoint records |
| endpoint_status | fail-closed adjudication: OBSERVED_STRICT_RECIST, OBSERVED_PCR, OBSERVED_OTHER_ENDPOINT, MULTIPLE_NATIVE_VALUES_REVIEW, ENDPOINT_DEFINITION_UNRESOLVED, LABEL_NOT_FOUND, MAPPING_NOT_FOUND, PUBLIC_FILE_ACCESS_CONFLICT |
| endpoint_verification_status | how endpoint rows joined to the patient (EXACT_SAMPLE:n / EXACT_PATIENT:n / NOT_LINKED) |

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
| clinical_endpoint_available | endpoint_status starts with OBSERVED_ |
| strict_auo_eligible | paired_valid AND endpoint_status == OBSERVED_STRICT_RECIST (PR/CR-vs-PD); the frozen manuscript A/U/O endpoint is the GSE91061 27-patient subset |

## Release files
`patients.parquet` (resource-level patient entries), `samples.parquet`,
`endpoints.parquet` (linked endpoint rows incl. recovery sources),
`provenance.parquet`, `duplicate_groups.csv`, `LABEL_RECOVERY_LEDGER.csv`,
`DATASET_SUMMARY.json`.
