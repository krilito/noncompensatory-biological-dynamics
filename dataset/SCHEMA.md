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

`gene_space/` (v0.2 canonical gene space) is documented in `gene_space/README.md`: its
feature maps, source metrics, gene-level union/core table, source overlap matrix,
cohort inventory and per-interval `PAIR_GENE_SPACE.csv.gz`. It adds no columns to the
MASTER schema and does not change any v0.1.x semantic field.
