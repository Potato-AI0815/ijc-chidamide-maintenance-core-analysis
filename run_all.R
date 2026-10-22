# =============================================================================
# IJC public core analysis entry point
#
# Scope: reconstructs the 47 treated + 23 platform-compatible comparator cohort
# from three coded frozen workbooks, then refits the propensity-score/overlap-
# weight pipeline, the time-dependent Cox models A/B/C, the existing sensitivity
# analyses (leave-one-comparator-out, common-support trimming, CR-only
# restriction) and the descriptive KM/CIF estimates.
#
# NOT in scope: figure generation, Table 1 display-row standardised mean
# differences, the final maintenance-period TEAE table (Supplementary Table S11,
# which incorporates six later author-confirmed hemorrhagic cystitis events that
# are absent from the frozen AE workbook), and the cytogenetic summary.
#
# Run from anywhere: Rscript path/to/run_all.R
# This script reads only the three coded frozen workbooks in data/.
# It reconstructs the current 47 + 23 cohort; no unprocessed patient records are
# needed.
# =============================================================================
invisible(suppressWarnings(Sys.setlocale('LC_CTYPE','.UTF-8')))
options(stringsAsFactors=FALSE,encoding='UTF-8')
required <- c('readxl','survival','cmprsk')
missing <- required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing))stop('Install required R packages first: ',paste(missing,collapse=', '))
args<-commandArgs(FALSE)
script<-args[grepl('^--file=',args)]
ROOT<-if(length(script))dirname(normalizePath(sub('^--file=','',script[1]),winslash='/'))else normalizePath('.')
OUT<-file.path(ROOT,'outputs');dir.create(OUT,showWarnings=FALSE)
source(file.path(ROOT,'code/model_functions_v2_9_2.R'),encoding='UTF-8')
source(file.path(ROOT,'code/absolute_rate_functions.R'),encoding='UTF-8')
source(file.path(ROOT,'code/sensitivity_functions.R'),encoding='UTF-8')
DAYS_PER_MONTH<-30.4375
read_sheet<-function(file,sheet)as.data.frame(readxl::read_excel(file.path(ROOT,'data',file),sheet=sheet,guess_max=10000))
CLINICAL<-'Chidamide_141_foundation_master_v11_final_frozen_deidentified.xlsx'
SAFETY<-'Chidamide_Safety_AE_FINAL_LOCKED_v1.1_deidentified.xlsx'
REFERENCE<-'STATISTICAL_FREEZE_RESULTS_v11_deidentified.xlsx'
pm<-read_sheet(CLINICAL,'Patient_Master_141')
endpoints<-read_sheet(CLINICAL,'Outcomes_141')
exposure<-read_sheet(CLINICAL,'Treatment_Exposure_47')
stopifnot(nrow(pm)==141,!anyDuplicated(pm$Patient_Key),nrow(endpoints)==141,!anyDuplicated(endpoints$Patient_Key))
pm$Treated_47<-as.integer(pm$Treatment_Group==1)
pm$Platform_Compatible_Control_23<-as.integer(pm$Treatment_Group==0 & pm$Calendar_Overlap_Flag==1 &
  pm$ATG_in_Conditioning==0 & pm$PTCy_CTX_in_Prophylaxis==1 & pm$FAB_Backbone_Indicator==1)
pm$Primary_Comparative_70<-as.integer(pm$Treated_47==1 | pm$Platform_Compatible_Control_23==1)
pm$Start_Delay_Days_Calculated<-exposure$Start_Delay_Days_Calculated[match(pm$Patient_Key,exposure$Patient_Key)]
d<-merge(pm,endpoints,by='Patient_Key',all.x=TRUE,sort=FALSE,suffixes=c('','_Outcome'))
d<-d[match(pm$Patient_Key,d$Patient_Key),]
stopifnot(sum(d$Treated_47)==47,sum(d$Platform_Compatible_Control_23)==23,sum(d$Primary_Comparative_70)==70)
write_utf8_csv(data.frame(Population=c('Source','Treated','Historical candidates','Calendar-overlap comparators','Platform-compatible comparators','Primary comparative'),
 N=c(nrow(d),sum(d$Treated_47),sum(d$Treatment_Group==0),sum(d$Treatment_Group==0&d$Calendar_Overlap_Flag==1),sum(d$Platform_Compatible_Control_23),sum(d$Primary_Comparative_70))),file.path(OUT,'cohort_counts.csv'))

# Same endpoint priority as v2.9.2: main endpoints use frozen month fields.
# For sensitivity GRFS, exact-date differences (now relative days) take priority
# over month fields. De-identification does not reconcile existing date/month discrepancies.
event_time<-function(flag,relative_days,months) {
  t<-rep(NA_real_,length(flag)); use<-!is.na(flag)&flag==1
  t[use]<-ifelse(is.finite(relative_days[use]),relative_days[use]/DAYS_PER_MONTH,months[use]);t
}
comp<-cbind(aGVHD=event_time(d$aGVHD_III_IV_Event,d$aGVHD_Day,rep(NA_real_,nrow(d))),
 cGVHD=event_time(d$Moderate_Severe_cGVHD_Event,rep(NA_real_,nrow(d)),d$Moderate_Severe_cGVHD_Months),
 Relapse=event_time(d$Relapse_Event,d$Relapse_Date_Derived_Days_From_HSCT,d$Relapse_Months),
 Death=event_time(d$OS_Event,d$OS_Date_Derived_Days_From_HSCT,d$OS_Months))
comp[!is.finite(comp)]<-NA_real_
d$GRFS_Event_Main<-as.integer(d$GRFS_Event);d$GRFS_Months_Main<-as_num(d$GRFS_Months)
d$GRFS_Event_Sensitivity<-as.integer(rowSums(!is.na(comp))>0)
d$GRFS_Months_Sensitivity<-apply(comp,1,function(x)if(all(is.na(x)))NA_real_ else min(x,na.rm=TRUE))
d$GRFS_Months_Sensitivity[d$GRFS_Event_Sensitivity==0]<-d$OS_Months[d$GRFS_Event_Sensitivity==0]
ctx<-new.env(parent=emptyenv());ctx$primary70<-d[d$Primary_Comparative_70==1,]
ctx$paths<-list(output_data=OUT,pipeline_logs=OUT);ctx$ps_lambda<-0.01;ctx$days_per_month<-DAYS_PER_MONTH
ctx$validation<-list();ctx$warnings<-character();set.seed(20260630)
invisible(pipeline_fit_propensity_score(ctx))
invisible(pipeline_build_balance_outputs(ctx))
invisible(pipeline_fit_survival_models(ctx))
write_utf8_csv(ctx$primary70,file.path(OUT,'primary_cohort_analysis.csv'))

# Absolute estimates: unweighted descriptive KM and CIF, not comparative effects.
d<-ctx$primary70
defs<-list(OS=c('OS_Months','OS_Event'),LFS=c('LFS_Months','LFS_Event'),
 'GRFS main definition'=c('GRFS_Months_Main','GRFS_Event_Main'),
 'GRFS sensitivity definition'=c('GRFS_Months_Sensitivity','GRFS_Event_Sensitivity'))
rows<-list()
for(g in c(1,0)) {
 z<-d[d$Treatment_Group==g,];gn<-if(g==1)'Treated' else 'Platform-compatible comparator'
 for(ep in names(defs))rows[[length(rows)+1]]<-km_points(z,defs[[ep]][1],defs[[ep]][2],c(12,24),ep,gn,'Frozen endpoint fields')
 cir<-as.integer(z$Competing_Risk_Code)
 rows[[length(rows)+1]]<-cif_points(z$CIR_NRM_Time_Months,ifelse(cir==1,1L,ifelse(cir==2,2L,0L)),c(12,24),'Relapse',gn,'months','Frozen event coding')
 rows[[length(rows)+1]]<-cif_points(z$CIR_NRM_Time_Months,ifelse(cir==2,1L,ifelse(cir==1,2L,0L)),c(12,24),'NRM',gn,'months','Frozen event coding')
 death<-ifelse(z$OS_Event==1,z$OS_Months,Inf);censor<-z$OS_Months
 ag<-make_competing_status(ifelse(z$aGVHD_III_IV_Event==1,z$aGVHD_Day,Inf),z$aGVHD_III_IV_Event,cbind(death*DAYS_PER_MONTH),censor*DAYS_PER_MONTH)
 rows[[length(rows)+1]]<-cif_points(ag$time,ag$status,100,'Grade III-IV aGVHD',gn,'days','Death as competing event')
 rel<-ifelse(z$Relapse_Event==1,z$CIR_NRM_Time_Months,Inf)
 cg<-make_competing_status(z$Moderate_Severe_cGVHD_Months,z$Moderate_Severe_cGVHD_Event,cbind(rel,death),censor)
 rows[[length(rows)+1]]<-cif_points(cg$time,cg$status,c(12,24),'Moderate/severe cGVHD',gn,'months','Direct CIF from HSCT; relapse and death competing')
}
write_utf8_csv(do.call(rbind,rows),file.path(OUT,'absolute_rates.csv'))
td<-ctx$td_intervals;ps<-ctx$ps_dataset;loo<-list();alt<-list()
for(ep in c('OS','LFS','GRFS_Main','GRFS_Sensitivity'))for(key in ps$Patient_Key[ps$Treatment_Group==0]) {
 cp<-td[td$Endpoint==ep&td$Patient_Key!=key,]
 ans<-safe_fit(fit_weighted_td(cp,ep,'Leave-one-comparator-out overlap-weighted time-dependent Cox','Original weights retained'))
 if(inherits(ans$value,'error'))stop(conditionMessage(ans$value))
 rr<-ans$value;rr$Omitted_Patient_Key<-key;rr$Warning<-ans$warning;loo[[length(loo)+1]]<-rr
}
write_utf8_csv(do.call(rbind,loo),file.path(OUT,'leave_one_comparator_out.csv'))
lo<-max(tapply(ps$Propensity_Score,ps$Treatment_Group,min));hi<-min(tapply(ps$Propensity_Score,ps$Treatment_Group,max))
specs<-list('PS common-support trimmed overlap-weighted time-dependent Cox'=ps$Patient_Key[ps$Propensity_Score>=lo&ps$Propensity_Score<=hi],
 'CR-only disease-risk restricted overlap-weighted time-dependent Cox'=ps$Patient_Key[ps$PreTx_Status=='CR'])
for(nm in names(specs))for(ep in c('OS','LFS','GRFS_Main','GRFS_Sensitivity')) {
 cp<-td[td$Endpoint==ep&td$Patient_Key%in%specs[[nm]],]
 ans<-safe_fit(fit_weighted_td(cp,ep,nm,'Original overlap weights retained; no PS refit within restriction'))
 if(inherits(ans$value,'error'))stop(conditionMessage(ans$value))
 rr<-ans$value;rr$Warning<-ans$warning;alt[[length(alt)+1]]<-rr
}
write_utf8_csv(do.call(rbind,alt),file.path(OUT,'alternative_adjustment.csv'))

# -----------------------------------------------------------------------------
# Frozen safety source summary -- NOT the final maintenance-period TEAE table.
#
# The published Supplementary Table S11 reports nine haemorrhagic cystitis cases
# (any grade). The frozen AE workbook below contains only the three cases present
# in the original AE listing; six further cases were later confirmed by the
# investigators from the transplant-complication record and are not part of this
# workbook. Running this section therefore yields HC = 3 and must never be
# presented as a reproduction of the published S11.
#
# Disabled by default so that a default run cannot be mistaken for the final
# safety table. Enable explicitly with the environment variable
# IJC_INCLUDE_FROZEN_SAFETY=true and label any output accordingly.
#
# Counting only: no new grading, attribution, SAE or DLT adjudication.
# -----------------------------------------------------------------------------
INCLUDE_FROZEN_SAFETY <- identical(tolower(Sys.getenv('IJC_INCLUDE_FROZEN_SAFETY')),'true')
if(INCLUDE_FROZEN_SAFETY) {
 ae<-read_sheet(SAFETY,'AE_Long_FINAL')
 stopifnot(length(unique(ae$Patient_Key))==47)
 count<-function(x)length(unique(x$Patient_Key))
 severe<-function(x)as.character(x$CTCAE_Grade)%in%c('3','4','3-4')|x$Grade_Group=='Grade 3-4'
 sr<-lapply(sort(unique(ae$AE_Term)),function(t) {
  z<-ae[ae$AE_Term==t,];r<-z[z$Related_to_Chidamide=='Yes',]
  data.frame(Term=t,Any_TEAE=count(z),Grade3=count(z[as.character(z$CTCAE_Grade)=='3',]),
  Grade4=count(z[as.character(z$CTCAE_Grade)=='4',]),Grade3_4_TEAE=count(z[severe(z),]),
  Any_TRAE=count(r),Grade3_4_TRAE=count(r[severe(r),]))
 })
 write_utf8_csv(do.call(rbind,sr),file.path(OUT,'frozen_safety_source_TERM_SUMMARY_NOT_FINAL_S11.csv'))
 write_utf8_csv(data.frame(Any_TEAE=count(ae),Any_TRAE=count(ae[ae$Related_to_Chidamide=='Yes',]),
  Severe_TEAE=count(ae[severe(ae),]),Severe_TRAE=count(ae[severe(ae)&ae$Related_to_Chidamide=='Yes',])),
  file.path(OUT,'frozen_safety_source_OVERVIEW_NOT_FINAL_S11.csv'))
}

hc<-list()
for(g in c(1,0))for(ep in c('CMV','HC')) {
 z<-d[d$Treatment_Group==g,]
 x<-if(ep=='CMV')z$CMV_Event else ifelse(z$Hemorrhagic_Cystitis_Event==0,0,ifelse(is.na(z$Hemorrhagic_Cystitis_Day),NA,as.integer(z$Hemorrhagic_Cystitis_Day<=100)))
 n<-sum(!is.na(x));k<-sum(x==1,na.rm=TRUE);ci<-binom.test(k,n)$conf.int
 hc[[length(hc)+1]]<-data.frame(Endpoint=ep,Group=if(g==1)'Treated' else 'Comparator',N=n,Events=k,Percent=100*k/n,Lower_95=100*ci[1],Upper_95=100*ci[2],Missing=sum(is.na(x)))
}
write_utf8_csv(do.call(rbind,hc),file.path(OUT,'frozen_complications.csv'))
# Treatment-delivery counts supported by structured fields only. No reconstruction
# of dose chronology from removed narrative records, and no new date adjudication.
dl<-read_sheet(CLINICAL,'Dose_Event_Log_47')
delivery<-data.frame(Metric=c('Patients reaching 20 mg per dose in the structured timeline',
 'Patients with structured dose-reduction events','Permanent discontinuation','Completed 365 days'),
 n=c(length(unique(dl$Patient_Key[!is.na(dl$Dose_mg_per_dose)&dl$Dose_mg_per_dose==20])),
 length(unique(dl$Patient_Key[dl$Event_Type=='Dose reduction'])),
 sum(exposure$Permanent_Discontinuation==1,na.rm=TRUE),sum(exposure$Completed_365_Days==1,na.rm=TRUE)),N=47)
write_utf8_csv(delivery,file.path(OUT,'structured_treatment_delivery.csv'))
# Reference outputs are validation targets ONLY; never fitted-data inputs.
ref<-read_sheet(REFERENCE,'TD_Cox');actual<-ctx$td_results[ctx$td_results$Term=='Exposure',]
idx<-match(paste(actual$Endpoint,actual$Model),paste(ref$Endpoint,ref$Model))
stopifnot(!anyNA(idx))
checks<-do.call(rbind,lapply(c('HR','CI_Lower_95','CI_Upper_95','P_Value','Events','N_Patients'),function(nm) {
 delta<-abs(actual[[nm]]-ref[[nm]][idx]);data.frame(Endpoint=actual$Endpoint,Model=actual$Model,
 Metric=nm,Recomputed=actual[[nm]],Reference=ref[[nm]][idx],Abs_Difference=delta,Pass=delta<1e-6)
}))
write_utf8_csv(checks,file.path(OUT,'reference_checks.csv'))
if(!all(checks$Pass))stop('Reference reproduction failed; see outputs/reference_checks.csv')
capture.output(sessionInfo(),file=file.path(OUT,'R_sessionInfo.txt'))
cat('PASS: cohort 47 + 23; all 12 main exposure-model estimates reproduced.\n')
cat('PASS: descriptive KM/CIF, leave-one-comparator-out and restriction analyses written.\n')
if(INCLUDE_FROZEN_SAFETY) {
  cat('WROTE: frozen safety source summary (HC=3). This is NOT the published Table S11.\n')
} else {
  cat('SKIPPED: frozen safety source summary (default). Set IJC_INCLUDE_FROZEN_SAFETY=true to enable.\n')
  cat('        The published Table S11 is not reproducible from this repository alone.\n')
}
