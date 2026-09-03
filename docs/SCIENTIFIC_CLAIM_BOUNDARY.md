# Scientific claim boundary

This file is the public claim ceiling for the accompanying manuscript:

> What paired biological measurements identify about biological change.

MELD-ICB is the frozen four-axis ruler used as a real-data evaluation object.
The repository does not claim a complete biological world model, a clinical
predictor, or a causal intervention effect.

Statuses are non-compensatory. A supported component cannot rescue an
unsupported, descriptive, conflict, or not-tested component.

| Component | Question | Status |
|-----------|----------|--------|
| B0 | Is the measurement system fixed? | PARTIAL / CONFLICT |
| B1 core | Does state transfer in PRJEB23709 and MORRISON-1-public? | SUPPORTED |
| B1 MGH biopsy | Does state transfer at MGH biopsy level? | NOT SUPPORTED |
| B1 MGH patient | Earliest eligible on-treatment biopsy per patient | SUPPORTIVE (different unit) |
| B2 | Does the state move? | SUPPORTED |
| B3 | Is movement outcome-specific? | NOT SUPPORTED |
| B4 | Does change add incremental information? | NOT SUPPORTED |
| B5 | Is biological meaning grounded? | DESCRIPTIVE |
| B6 | Where does interpretation stop? | ENVELOPE DECLARED |
| B-INT | Is movement intervention-specific? | NOT TESTED |

## B0

The public implementation has one canonical preprocessing path:
`PRE_MEAN_THEN_AXIS_Z`, `ddof=1`, duplicate symbols collapsed by mean,
coverage fail-closed. Coefficients and the threshold are frozen:

`-0.676 T + 0.173 E + 0 X + 1.135 C - 0.666`, threshold `0.735`.

The manuscript also reports a historical preprocessing lineage conflict.
That conflict is a historical-path fact, not a second live preprocessing
switch in this repository. `source_data/figure2_auc.csv` retains both
`auc_canonical` / `p_canonical` and `auc_historical` / `p_historical`.
Do not delete the historical columns to make B0 look clean.

## Forbidden upgrades

- Do not invert an AUC because it is below 0.5.
- Do not write MGH biopsy-level transfer as supported.
- Do not write B3 or B4 as supported.
- Do not write B5 as patient-level inference or as a cross-substrate causal switch.
- Do not write B6 context/comparison/superseded cohorts as external validation.
- Do not redistribute IMvigor210 sample-level scores. Accession and the
  paper-reported summary AUC/CI/P are the public objects.
- Do not write B-INT as “no effect”. It was not tested.
- Do not treat pooled 85-biopsy AUC 0.771 as a cohort-level or meta-analytic result.
- Do not treat record-level cosine and cohort-mean vector cosine as one metric.
