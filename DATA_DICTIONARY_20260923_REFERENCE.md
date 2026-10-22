# Data dictionary and transformation rules

The original field names are retained wherever an allowed source field is copied unchanged. Blank cells remain missing; they are not replaced with zeros. `workbook_inventory.csv` lists every retained sheet and its dimensions.

## Clinical workbook

- `Patient_Master_141`: one row per coded research participant. `Treatment_Group` is 1 for the prospective maintenance cohort and 0 for historical comparator candidates. `Patient_Key` is the same pseudonymous join key across workbooks.
- Cohort selection: treated subjects are retained; comparators require `Calendar_Overlap_Flag = 1`, `ATG_in_Conditioning = 0`, `PTCy_CTX_in_Prophylaxis = 1`, and `FAB_Backbone_Indicator = 1`. Expected primary-cohort counts are 47 and 23. The calendar flag is retained without exact transplant dates.
- `FAB_Backbone_Indicator` is the numeric result of the existing pipeline's `grepl("FAB", Conditioning_Regimen_Raw, fixed = TRUE)` rule. It replaces that raw-text dependency; it does not mean all drugs or all care were identical. Both primary cohorts received chidamide-containing conditioning according to the current author confirmation. The comparison concerns post-transplant maintenance.
- `Transplant_Year` is retained for adjustment. Age, sex, diagnosis, pretreatment burden, remission/MRD classification, donor type, CD34 dose, prophylaxis and recorded eligibility fields are copied without reclassification. A retained eligibility flag is not a new clinical verification.
- `HLA_Match` is generalized to the exact dichotomy already used by the PS model and Table 1: `12/12` versus `Other/not recorded as 12/12`, preserving missing values. This removes a date-formatted source string and a donor-relationship label without changing any model/table classification. It does not clinically adjudicate those original entries or reconstruct detailed HLA matching.
- `Treatment_Exposure_47`: one row per treated participant. Start-delay fields are days from HSCT; `Starting_Dose_mg_per_dose` and `Maximum_Maintenance_Dose_mg` are mg per administration. Course-status flags remain as stored. Exposure in the time-dependent model switches on at `Start_Delay_Days_Calculated` and does not switch off after discontinuation.
- `Dose_Event_Log_47`: structured dose events only. `Day_Post_HSCT` is relative time; dose is mg per administration. Raw event narratives and calendar dates are excluded. Structured records and broad dose-reduction flags measure different things; they must not be silently equated.
- `Outcomes_141`: frozen clinical endpoints for all141 source subjects. `_Event` fields retain their source binary coding. `Competing_Risk_Code` is 0 for censoring, 1 for relapse and 2 for NRM. `_Months` fields retain their frozen values. Acute-GVHD, CMV, cystitis and infection `_Day` fields are days from HSCT.
- `OS_Date_Derived_Days_From_HSCT` and `Relapse_Date_Derived_Days_From_HSCT` replace the actual-date inputs with `event date - transplant date`, using the same Excel date origin as the existing pipeline. No actual calendar dates or offsets permitting their reconstruction are supplied.
- Derived main GRFS uses frozen `GRFS_Months`/`GRFS_Event`. Derived sensitivity GRFS is the first qualifying grade III–IV acute-GVHD, moderate/severe chronic-GVHD, relapse or death component under the existing time-priority rule. Date-derived relative days take priority over frozen months for relapse/death, as in the current code. The days/month conversion constant is 30.4375 where that code uses it.
- Descriptive moderate/severe cGVHD is direct CIF from HSCT with relapse/death as competing events; it is not a day-100 landmark estimate.

Excluded clinical sheets: duplicate `Patient_Master_124` and stale `Outcomes_124`; narrative-rich `Treatment_Baseline_Detail_47`; `Source_Treatment_47`, `Source_Control_94`; redundant eligibility sheet; and QC, decision, correction, deleted-row and original README logs. Their removal is documented rather than silently retaining inaccessible hidden data.

## Safety workbook

- `Subject_Level_FINAL`:47 coded participants and frozen summary flags/counts.
- `AE_Long_FINAL`: existing patient–event-term analysis rows. `CTCAE_Grade` and `Grade_Group` are retained as recorded/derived in the freeze, including combined grade1–2 categories. No new grades are assigned.
- `Related_to_Chidamide = Yes` selects investigator-attributed TRAE rows under the frozen rules. TEAE and TRAE summaries are counted separately, once per participant per term. Grade3/4 combined flags must not be treated as a formal SAE or DLT classification.
- `Numeric_Value_Used`, `Numeric_Unit` and `Grade_Method` preserve the coded grading basis; they do not include the removed raw clinical narrative. `Supplemented_AE` is retained as an existing provenance flag.
- `Safety_Overview`, `Table2_FINAL`, `Supplementary_AE_Table` and `Grading_Rules` retain their aggregate cell content/rules. These include the original frozen source classifications, not newly adjudicated values.
- Excluded: raw AE text, related-AE narrative, dose blood-count text, source notes and the patient-level supplemented/relatedness audit sheets. The source workbook's internal README is not copied.

## Statistical-results workbook

Retained: `Cohorts`, `TD_Cox`, `ESS`, `Balance`, `Table1_Source`, `Table2_Source`, `Absolute_KM`, `Competing_Risk`, `PS_Model`. Values are reference outputs, not substitutes for patient-level estimation. The internal `Validation` log is omitted. Retention of a historical result sheet does not promote every analysis in it to a current primary endpoint or change the current manuscript's descriptive interpretation of competing-risk outcomes.

## R specifications

- Primary-cohort median imputation and standardization follow existingv2.9.2 functions. Categorical levels are unchanged.
- Logistic PS stability fallback uses ridge lambda0.01, with an unpenalized intercept. Weights are1−PS for treated subjects andPS for comparators.
- Cox Models A/B/C retain existing covariates, Breslow ties and robust patient-clustered variance. Model B is the main overlap-weighted time-dependent model.
- Landmark days30/60/100 fix exposure at the landmark. Leave-one-comparator-out, empirical common-support trimming and CR-only analyses retain original primary-cohort weights, without refitting PS within each subset.
- KM confidence intervals use the existing log–log method. CIF intervals use the existing normal approximation clipped to[0,1]. CMV and day-100 cystitis summaries use the current freeze only, with exact binomial intervals; the comparator cystitis denominator is23, not20.

These are descriptions of the supplied code, not additional methodological claims or new analyses.
