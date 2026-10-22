# Core analysis code — early chidamide maintenance after allogeneic HSCT

Analysis code supporting the manuscript:

> **Early chidamide maintenance after allogeneic hematopoietic stem cell transplantation in high-risk myeloid malignancies: a prospective cohort with time-dependent external-control analysis**

This repository contains the statistical analysis code for the primary comparative
analysis. It does **not** contain patient-level data.

---

## What this code does

`run_all.R` reconstructs the comparative cohort and refits every model reported in
the main text and the sensitivity analyses of Supplementary Tables S4A and S4B.

In scope:

- Reconstruction of the **47 treated + 23 platform-compatible comparator** cohort
  (primary comparative cohort, n = 70) from a 141-patient source cohort.
- Propensity-score model (ridge-stabilised logistic regression, fixed penalty 0.01)
  and overlap weights (`treated = 1 − PS`, `comparator = PS`).
- Effective sample size and balance summaries.
- Time-dependent Cox models **A** (unweighted, adjusted), **B** (overlap-weighted,
  exposure only — the main exploratory comparison) and **C** (overlap-weighted,
  Model A covariates), with Breslow ties and patient-clustered robust variance.
- Existing sensitivity analyses: leave-one-comparator-out (23 fits per endpoint),
  common-support trimming, and CR-only restriction. Overlap weights are retained
  and the propensity score is not refitted within restrictions.
- Descriptive, unweighted Kaplan–Meier and cumulative-incidence estimates at
  12 and 24 months (and day 100 for grade III–IV aGVHD). These are descriptive and
  are not adjusted treatment effects.

Not in scope — do not expect these files to reproduce the following:

| Not reproduced | Reason |
|---|---|
| **Supplementary Table S11** (maintenance-period TEAE) | The published table reports nine haemorrhagic cystitis cases (any grade). The frozen AE workbook used here contains only the three cases from the original AE listing; six further cases were confirmed later by the investigators from the transplant-complication record. Those additions, and the case-level decisions behind them, are held in controlled data and are deliberately not part of this repository. |
| Table 1 display-row standardised mean differences | Display-table row-level SMDs are not identical to model-matrix term SMDs. |
| Table 2 / Supplementary Table S10 (treatment-related AE) | Derived from the frozen safety workbook through a separate attribution workflow. |
| Supplementary Tables S5 (cytogenetics), S6 (treatment delivery), S9 (NRM causes) | Separate input sources and/or investigator-confirmed conventions. |
| All figures (1–4, S1–S4) | Figure production uses a separate, multi-stage plotting chain that is not part of this release. |

A default run therefore **does not** write any safety-summary output, so that a
default run cannot be mistaken for a reproduction of the published TEAE table.

---

## Requirements

- **R 4.6.1** (verified). Other recent R versions are expected to work but were not tested.
- R packages: `readxl`, `survival`, `cmprsk`.

The versions below are the ones actually used for the verification run recorded in
`validation/`. They are this project's *reproduction* environment, not the original
historical analysis environment.

| Package | Version used |
|---|---|
| survival | 3.8.6 |
| cmprsk | 2.2.12 |
| readxl | 1.6.0 |

A `sessionInfo()` capture is written to `outputs/R_sessionInfo.txt` on every run.

If you want a reproducible package library, initialise `renv` and commit the
resulting `renv.lock`. Note that `renv` pins package versions only; it does not
capture the R version itself or system-level dependencies, so record those
separately.

---

## Input data

Input data are **not** included in this repository and are **not** publicly
available.

`run_all.R` expects three coded (de-identified) frozen workbooks in a `data/`
directory:

| File | Contents |
|---|---|
| `Chidamide_141_foundation_master_v11_final_frozen_deidentified.xlsx` | Clinical source; sheets `Patient_Master_141`, `Outcomes_141`, `Treatment_Exposure_47`, `Dose_Event_Log_47` |
| `Chidamide_Safety_AE_FINAL_LOCKED_v1.1_deidentified.xlsx` | Frozen AE listing; sheet `AE_Long_FINAL` |
| `STATISTICAL_FREEZE_RESULTS_v11_deidentified.xlsx` | Frozen reference estimates used as validation targets only; sheet `TD_Cox` |

Required columns and sheet schemas are documented in
`DATA_DICTIONARY_20260923_REFERENCE.md`. Please read it before requesting data —
it describes exactly what the code needs.

Individual-level patient data are protected by ethical and institutional
regulations. Access requires institutional approval and is not guaranteed. The
corresponding author can be contacted for further information.

The workbooks remain coded patient-level data even after de-identification. They
are not aggregate data and must not be redistributed.

---

## Running

```sh
# place the three authorised workbooks in ./data/
Rscript run_all.R
```

Outputs are written to `./outputs/`:

| Output | Content |
|---|---|
| `cohort_counts.csv` | Source 141 / treated 47 / historical 94 / calendar-overlap 77 / platform-compatible 23 / primary 70 |
| `primary_cohort_analysis.csv` | Patient-level analysis dataset — **contains coded patient records; keep private** |
| `absolute_rates.csv` | Descriptive KM/CIF estimates |
| `leave_one_comparator_out.csv` | 23 leave-one-out fits × 4 endpoints |
| `alternative_adjustment.csv` | Common-support trimming and CR-only restriction |
| `reference_checks.csv` | 72 recomputed-versus-frozen comparisons; the run aborts if any differ by ≥ 1e-6 |
| `frozen_complications.csv` | Descriptive CMV reactivation and grade 2–4 haemorrhagic cystitis by day 100 |
| `structured_treatment_delivery.csv` | Treatment-delivery counts from structured fields only |
| `R_sessionInfo.txt` | Environment capture |

`outputs/primary_cohort_analysis.csv` contains coded patient-level rows. **Do not
run this script into a publicly synchronised directory.**

### Optional: frozen safety source summary

```sh
IJC_INCLUDE_FROZEN_SAFETY=true Rscript run_all.R
```

This additionally writes `frozen_safety_source_TERM_SUMMARY_NOT_FINAL_S11.csv` and
`frozen_safety_source_OVERVIEW_NOT_FINAL_S11.csv`. These describe the frozen AE
workbook only and give haemorrhagic cystitis any-grade **3**, not the published 9.
They are named explicitly so they cannot be mistaken for Supplementary Table S11.

---

## Verification

`validation/` holds reference outputs from the verification run, containing
aggregate results only (no patient-level rows):

- `reference_checks.csv` — 72 rows, all passing, maximum absolute difference
  ≈ 4.2 × 10⁻¹⁵.
- `cohort_counts.csv` — cohort reconstruction.
- `absolute_rates.csv` — descriptive estimates.

`run_all.R` re-derives `reference_checks.csv` on every run and calls `stop()` if any
comparison fails, so a successful run is itself the check.

---

## License

See `LICENSE`.

## Contact

Corresponding authors: Yuhuan Zheng (zhengyuhuan@scu.edu.cn) and Jie Ji
(jieji@scu.edu.cn), Department of Hematology / Institute of Hematology, West China
Hospital, Sichuan University, Chengdu, China.
