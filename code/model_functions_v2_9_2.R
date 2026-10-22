# Function-only extraction of the existing v2.9.2 model/balance code. Statistical formulas unchanged.

add_validation <- 
function (ctx, check_id, category, item, expected, observed, status = "PASS", details = "") 
{
    ctx$validation[[length(ctx$validation) + 1L]] <- data.frame(Check_ID = check_id, Category = category, Item = item, Expected = as.character(expected), 
        Observed = as.character(observed), Status = status, Details = as.character(details), stringsAsFactors = FALSE)
    invisible(ctx)
}

append_log <- 
function (ctx, text) 
{
    path <- file.path(ctx$paths$pipeline_logs, "REPRODUCIBLE_PIPELINE_LOG.txt")
    line <- paste0(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), " | ", text, "\n")
    cat(line, file = path, append = TRUE)
    invisible(path)
}

as_num <- 
function (x) 
suppressWarnings(as.numeric(x))

binary_summary <- 
function (ind, treat, w) 
{
    data.frame(Treated_Before = paste0(sum(ind[treat == 1]), " (", fmt_num(100 * mean(ind[treat == 1]), 1), "%)"), Control_Before = paste0(sum(ind[treat == 
        0]), " (", fmt_num(100 * mean(ind[treat == 0]), 1), "%)"), Treated_After = paste0(fmt_num(100 * weighted_mean_safe(ind[treat == 
        1], w[treat == 1]), 1), "%"), Control_After = paste0(fmt_num(100 * weighted_mean_safe(ind[treat == 0], w[treat == 
        0]), 1), "%"), stringsAsFactors = FALSE)
}

build_counting_process <- 
function (d, time_var, event_var, dpm) 
{
    rows <- vector("list", nrow(d) * 2L)
    k <- 0L
    for (i in seq_len(nrow(d))) {
        end <- as_num(d[[time_var]][i])
        ev <- as.integer(d[[event_var]][i] == 1)
        if (!is.finite(end) || end <= 0) 
            next
        treated <- d$Treatment_Group[i] == 1
        start <- if (treated) 
            as_num(d$Start_Delay_Days_Calculated[i])/dpm
        else Inf
        base <- data.frame(Patient_Key = d$Patient_Key[i], Treatment_Group = d$Treatment_Group[i], Propensity_Score = d$Propensity_Score[i], 
            Overlap_Weight = d$Overlap_Weight[i], Age10 = (as_num(d$Age[i]) - mean(as_num(d$Age), na.rm = TRUE))/10, NonCR = as.integer(d$PreTx_Status[i] != 
                "CR"), MRD_Risk = as.integer(d$MRD_Status[i] != "Negative"), HID = as.integer(d$Donor_Type[i] == "HID"), 
            TxYearCentered = as_num(d$Transplant_Year[i]) - mean(as_num(d$Transplant_Year), na.rm = TRUE), stringsAsFactors = FALSE)
        if (treated && is.finite(start) && start > 0 && start < end) {
            k <- k + 1L
            rows[[k]] <- cbind(base, tstart = 0, tstop = start, event = 0L, Exposure = 0L)
            k <- k + 1L
            rows[[k]] <- cbind(base, tstart = start, tstop = end, event = ev, Exposure = 1L)
        }
        else {
            exposure <- as.integer(treated && is.finite(start) && start <= 0)
            k <- k + 1L
            rows[[k]] <- cbind(base, tstart = 0, tstop = end, event = ev, Exposure = exposure)
        }
    }
    do.call(rbind, rows[seq_len(k)])
}

extract_cox <- 
function (fit, endpoint, model_name, time_origin, exposure_definition, n_patients = NULL) 
{
    s <- summary(fit)
    co <- as.data.frame(s$coefficients)
    if (!nrow(co)) 
        return(data.frame())
    se_col <- if ("robust se" %in% names(co)) 
        "robust se"
    else "se(coef)"
    p_col <- grep("^Pr\\(", names(co), value = TRUE)[1]
    out <- data.frame(Endpoint = endpoint, Model = model_name, Term = rownames(co), Log_HR = co$coef, Robust_SE = co[[se_col]], 
        HR = exp(co$coef), CI_Lower_95 = exp(co$coef - 1.96 * co[[se_col]]), CI_Upper_95 = exp(co$coef + 1.96 * co[[se_col]]), 
        P_Value = co[[p_col]], N_Patients = if (is.null(n_patients)) 
            fit$n
        else n_patients, N_Intervals = fit$n, Events = fit$nevent, Time_Origin = time_origin, Exposure_Definition = exposure_definition, 
        Variance = "Robust sandwich clustered by Patient_Key", stringsAsFactors = FALSE)
    rownames(out) <- NULL
    out
}

fit_landmark <- 
function (d, time_var, event_var, endpoint, landmark_day, dpm) 
{
    landmark_month <- landmark_day/dpm
    keep <- is.finite(as_num(d[[time_var]])) & as_num(d[[time_var]]) > landmark_month
    x <- d[keep, , drop = FALSE]
    if (nrow(x) < 10L) 
        return(data.frame())
    x$Time_From_Landmark <- as_num(x[[time_var]]) - landmark_month
    x$Event_After_Landmark <- as.integer(x[[event_var]] == 1)
    x$Started_By_Landmark <- as.integer(x$Treatment_Group == 1 & as_num(x$Start_Delay_Days_Calculated) <= landmark_day)
    fit <- survival::coxph(survival::Surv(Time_From_Landmark, Event_After_Landmark) ~ Started_By_Landmark + survival::cluster(Patient_Key), 
        data = x, weights = Overlap_Weight, robust = TRUE, ties = "breslow", model = TRUE, x = TRUE)
    out <- extract_cox(fit, endpoint, paste0("Overlap-weighted landmark day ", landmark_day), paste0("Day ", landmark_day, 
        " after HSCT"), paste0("Fixed status at landmark: maintenance started by day ", landmark_day), n_patients = nrow(x))
    out$Landmark_Day <- landmark_day
    out
}

fit_td_endpoint <- 
function (d, time_var, event_var, endpoint, dpm) 
{
    cp <- build_counting_process(d, time_var, event_var, dpm)
    fits <- list()
    fits$Adjusted_unweighted <- survival::coxph(survival::Surv(tstart, tstop, event) ~ Exposure + Age10 + NonCR + MRD_Risk + 
        HID + TxYearCentered + survival::cluster(Patient_Key), data = cp, ties = "breslow", robust = TRUE, model = TRUE, 
        x = TRUE)
    fits$Overlap_weighted <- survival::coxph(survival::Surv(tstart, tstop, event) ~ Exposure + survival::cluster(Patient_Key), 
        data = cp, weights = Overlap_Weight, ties = "breslow", robust = TRUE, model = TRUE, x = TRUE)
    fits$Overlap_weighted_residual_adjusted <- survival::coxph(survival::Surv(tstart, tstop, event) ~ Exposure + Age10 + 
        NonCR + MRD_Risk + HID + TxYearCentered + survival::cluster(Patient_Key), data = cp, weights = Overlap_Weight, ties = "breslow", 
        robust = TRUE, model = TRUE, x = TRUE)
    result <- do.call(rbind, lapply(names(fits), function(nm) {
        extract_cox(fits[[nm]], endpoint, nm, "Allogeneic HSCT", "Time-dependent maintenance exposure switches from 0 to 1 on the actual start date", 
            n_patients = length(unique(cp$Patient_Key)))
    }))
    list(intervals = cp, fits = fits, result = result)
}

fmt_num <- 
function (x, digits = 3L) 
{
    ifelse(is.finite(x), formatC(x, digits = digits, format = "f"), "")
}

impute_numeric_with_indicator <- 
function (x) 
{
    x <- as_num(x)
    med <- stats::median(x, na.rm = TRUE)
    if (!is.finite(med)) 
        med <- 0
    list(value = ifelse(is.na(x), med, x), missing = as.integer(is.na(x)), median = med)
}

median_iqr_text <- 
function (x) 
{
    x <- as_num(x)
    q <- stats::quantile(x, c(0.25, 0.5, 0.75), na.rm = TRUE, names = FALSE)
    paste0(fmt_num(q[2], 1), " [", fmt_num(q[1], 1), ", ", fmt_num(q[3], 1), "]")
}

pipeline_build_balance_outputs <- 
function (ctx) 
{
    d <- ctx$primary70
    w <- d$Overlap_Weight
    treat <- as.integer(d$Treatment_Group)
    rows <- list()
    continuous <- list(Age = as_num(d$Age), Diagnosis_to_HSCT_Days = as_num(d$Diagnosis_to_HSCT_Days_Calculated), PreTx_Treatment_Cycles = as_num(d$PreTx_Treatment_Cycles), 
        PreTx_BM_Blast_Pct = as_num(d$PreTx_BM_Blast_Pct), PreTx_FCM_Pct = as_num(d$PreTx_FCM_Pct), CD34_Dose_10e6_per_kg = as_num(d$CD34_Dose_10e6_per_kg))
    labels <- c(Age = "Age, years", Diagnosis_to_HSCT_Days = "Diagnosis-to-HSCT interval, days", PreTx_Treatment_Cycles = "Pretransplant treatment cycles", 
        PreTx_BM_Blast_Pct = "Pretransplant bone marrow blasts, %", PreTx_FCM_Pct = "Pretransplant flow-cytometric disease, %", 
        CD34_Dose_10e6_per_kg = "CD34+ cell dose, 10^6/kg")
    balance_term <- c(Age = "Age_z", Diagnosis_to_HSCT_Days = "Diagnosis_to_HSCT_Days_z", PreTx_Treatment_Cycles = "PreTx_Treatment_Cycles_z", 
        PreTx_BM_Blast_Pct = "PreTx_BM_Blast_Pct_z", PreTx_FCM_Pct = "PreTx_FCM_Pct_z", CD34_Dose_10e6_per_kg = "CD34_Dose_z")
    for (nm in names(continuous)) {
        x <- continuous[[nm]]
        imputed <- impute_numeric_with_indicator(x)$value
        before <- smd_one(imputed, treat)
        denom <- sqrt((stats::var(imputed[treat == 1]) + stats::var(imputed[treat == 0]))/2)
        after <- smd_one(imputed, treat, w, denom)
        rows[[length(rows) + 1L]] <- data.frame(Variable = labels[nm], Level = "", Variable_Type = "Continuous", Treated_Before = median_iqr_text(x[treat == 
            1]), Control_Before = median_iqr_text(x[treat == 0]), SMD_Before = before, Treated_After = weighted_mean_sd_text(imputed[treat == 
            1], w[treat == 1]), Control_After = weighted_mean_sd_text(imputed[treat == 0], w[treat == 0]), SMD_After = after, 
            PS_Covariate = TRUE, PS_Model_Term = balance_term[nm], Missing_Treated = sum(is.na(x[treat == 1])), Missing_Control = sum(is.na(x[treat == 
                0])), stringsAsFactors = FALSE)
    }
    cat_vars <- list(Sex = ifelse(is.na(d$Sex), "Missing", d$Sex), Diagnosis = ifelse(is.na(d$Diagnosis), "Missing", d$Diagnosis), 
        PreTx_Status = ifelse(is.na(d$PreTx_Status), "Missing", d$PreTx_Status), MRD_Status = ifelse(is.na(d$MRD_Status), 
            "Missing", d$MRD_Status), Donor_Type = ifelse(is.na(d$Donor_Type), "Missing", d$Donor_Type), HLA_12_12 = ifelse(d$HLA_Match == 
            "12/12", "12/12", "Other/not recorded as 12/12"), Transplant_Year = as.character(d$Transplant_Year), MMF_in_Prophylaxis = ifelse(as_num(d$MMF_in_Prophylaxis) == 
            1, "Yes", "No"))
    cat_labels <- c(Sex = "Sex", Diagnosis = "Diagnosis", PreTx_Status = "Pretransplant disease status", MRD_Status = "Pretransplant MRD status", 
        Donor_Type = "Donor type", HLA_12_12 = "HLA match", Transplant_Year = "Transplant year", MMF_in_Prophylaxis = "MMF in GVHD prophylaxis")
    for (nm in names(cat_vars)) {
        x <- as.character(cat_vars[[nm]])
        for (lev in sort(unique(x))) {
            ind <- as.integer(x == lev)
            bs <- binary_summary(ind, treat, w)
            before <- smd_one(ind, treat)
            denom <- sqrt((stats::var(ind[treat == 1]) + stats::var(ind[treat == 0]))/2)
            after <- smd_one(ind, treat, w, denom)
            rows[[length(rows) + 1L]] <- data.frame(Variable = cat_labels[nm], Level = lev, Variable_Type = "Categorical", 
                Treated_Before = bs$Treated_Before, Control_Before = bs$Control_Before, SMD_Before = before, Treated_After = bs$Treated_After, 
                Control_After = bs$Control_After, SMD_After = after, PS_Covariate = TRUE, PS_Model_Term = paste0(nm, ":", 
                  lev), Missing_Treated = 0L, Missing_Control = 0L, stringsAsFactors = FALSE)
        }
    }
    missing_vars <- list(Missing_PreTx_Treatment_Cycles = is.na(d$PreTx_Treatment_Cycles), Missing_PreTx_BM_Blast_Pct = is.na(d$PreTx_BM_Blast_Pct), 
        Missing_PreTx_FCM_Pct = is.na(d$PreTx_FCM_Pct), Missing_CD34_Dose = is.na(d$CD34_Dose_10e6_per_kg))
    for (nm in names(missing_vars)) {
        ind <- as.integer(missing_vars[[nm]])
        bs <- binary_summary(ind, treat, w)
        denom <- sqrt((stats::var(ind[treat == 1]) + stats::var(ind[treat == 0]))/2)
        rows[[length(rows) + 1L]] <- data.frame(Variable = paste0("Missingness indicator: ", sub("^Missing_", "", nm)), Level = "Missing", 
            Variable_Type = "Missingness indicator", Treated_Before = bs$Treated_Before, Control_Before = bs$Control_Before, 
            SMD_Before = smd_one(ind, treat), Treated_After = bs$Treated_After, Control_After = bs$Control_After, SMD_After = smd_one(ind, 
                treat, w, denom), PS_Covariate = TRUE, PS_Model_Term = nm, Missing_Treated = 0L, Missing_Control = 0L, stringsAsFactors = FALSE)
    }
    table1 <- do.call(rbind, rows)
    table1$Absolute_SMD_Before <- abs(table1$SMD_Before)
    table1$Absolute_SMD_After <- abs(table1$SMD_After)
    patient_cols <- c("Patient_Key", "Treatment_Group", "Propensity_Score", "Overlap_Weight", "Age", "Sex", "Diagnosis", 
        "Diagnosis_to_HSCT_Days_Calculated", "PreTx_Treatment_Cycles", "PreTx_BM_Blast_Pct", "PreTx_FCM_Pct", "PreTx_Status", 
        "MRD_Status", "Donor_Type", "HLA_Match", "Transplant_Year", "CD34_Dose_10e6_per_kg", "MMF_in_Prophylaxis")
    patient_cols <- intersect(patient_cols, names(d))
    patient_data <- d[, patient_cols, drop = FALSE]
    trace <- data.frame(Table1_Variable = c(labels, cat_labels, "Four numeric missingness indicators"), Locked_Source_Field = c("Age", 
        "Diagnosis_to_HSCT_Days_Calculated", "PreTx_Treatment_Cycles", "PreTx_BM_Blast_Pct", "PreTx_FCM_Pct", "CD34_Dose_10e6_per_kg", 
        "Sex", "Diagnosis", "PreTx_Status", "MRD_Status", "Donor_Type", "HLA_Match", "Transplant_Year", "MMF_in_Prophylaxis", 
        "Missingness derived before median imputation"), Population = "Primary comparative cohort (47 treated; 23 controls)", 
        Weighting = "Single overlap weight", stringsAsFactors = FALSE)
    ctx$table1_source <- table1
    ctx$table1_patient_data <- patient_data
    ctx$table1_traceability <- trace
    write_utf8_csv(table1, file.path(ctx$paths$output_data, "table1_source_data.csv"))
    write_utf8_csv(patient_data, file.path(ctx$paths$output_data, "table1_patient_level_data.csv"))
    write_utf8_csv(trace, file.path(ctx$paths$output_data, "table1_traceability.csv"))
    max_after <- max(table1$Absolute_SMD_After[is.finite(table1$Absolute_SMD_After)])
    add_validation(ctx, "RPL-T1-001", "Table 1", "Primary comparative cohort", "47 treated and 23 controls", paste0(sum(treat == 
        1), " treated; ", sum(treat == 0), " controls"), if (sum(treat == 1) == 47L && sum(treat == 0) == 23L) 
        "PASS"
    else "FAIL")
    add_validation(ctx, "RPL-T1-002", "Table 1", "Complete pre/post SMD", "No missing SMD for any PS covariate row", paste0(nrow(table1), 
        " rows; missing before=", sum(is.na(table1$SMD_Before)), "; missing after=", sum(is.na(table1$SMD_After))), if (!any(is.na(table1$SMD_Before)) && 
        !any(is.na(table1$SMD_After))) 
        "PASS"
    else "FAIL")
    add_validation(ctx, "RPL-T1-003", "Table 1", "Publication-table row maximum absolute SMD after weighting", "<=0.10", 
        fmt_num(max_after, 6), if (max_after <= 0.1) 
            "PASS"
        else if (max_after <= 0.2) 
            "WARN"
        else "FAIL")
    append_log(ctx, paste0("Built Table 1 source with complete unweighted/weighted SMDs for all PS covariates. ", "Max |SMD| after=", 
        fmt_num(max_after, 4), "."))
    ctx
}

pipeline_fit_propensity_score <- 
function (ctx) 
{
    d <- ctx$primary70
    y <- as.integer(d$Treatment_Group)
    age <- impute_numeric_with_indicator(d$Age)
    dx <- impute_numeric_with_indicator(d$Diagnosis_to_HSCT_Days_Calculated)
    cycles <- impute_numeric_with_indicator(d$PreTx_Treatment_Cycles)
    bm <- impute_numeric_with_indicator(d$PreTx_BM_Blast_Pct)
    fcm <- impute_numeric_with_indicator(d$PreTx_FCM_Pct)
    cd34 <- impute_numeric_with_indicator(d$CD34_Dose_10e6_per_kg)
    z_age <- standardize_fixed(age$value)
    z_dx <- standardize_fixed(dx$value)
    z_cycles <- standardize_fixed(cycles$value)
    z_bm <- standardize_fixed(bm$value)
    z_fcm <- standardize_fixed(fcm$value)
    z_cd34 <- standardize_fixed(cd34$value)
    psdf <- data.frame(Age_z = z_age$value, Sex = factor(ifelse(is.na(d$Sex), "Missing", d$Sex)), Diagnosis = factor(ifelse(is.na(d$Diagnosis), 
        "Missing", d$Diagnosis)), Diagnosis_to_HSCT_Days_z = z_dx$value, PreTx_Treatment_Cycles_z = z_cycles$value, PreTx_BM_Blast_Pct_z = z_bm$value, 
        PreTx_FCM_Pct_z = z_fcm$value, PreTx_Status = factor(ifelse(is.na(d$PreTx_Status), "Missing", d$PreTx_Status)), MRD_Status = factor(ifelse(is.na(d$MRD_Status), 
            "Missing", d$MRD_Status)), Donor_Type = factor(ifelse(is.na(d$Donor_Type), "Missing", d$Donor_Type)), HLA_12_12 = as.integer(d$HLA_Match == 
            "12/12"), Transplant_Year = factor(d$Transplant_Year), CD34_Dose_z = z_cd34$value, MMF_in_Prophylaxis = as.integer(d$MMF_in_Prophylaxis), 
        Missing_PreTx_Treatment_Cycles = cycles$missing, Missing_PreTx_BM_Blast_Pct = bm$missing, Missing_PreTx_FCM_Pct = fcm$missing, 
        Missing_CD34_Dose = cd34$missing, check.names = FALSE)
    formula_text <- paste("Treatment_Group ~ Age_z + Sex + Diagnosis + Diagnosis_to_HSCT_Days_z +", "PreTx_Treatment_Cycles_z + PreTx_BM_Blast_Pct_z + PreTx_FCM_Pct_z +", 
        "PreTx_Status + MRD_Status + Donor_Type + HLA_12_12 + Transplant_Year +", "CD34_Dose_z + MMF_in_Prophylaxis + four numeric missingness indicators")
    mm <- stats::model.matrix(~., data = psdf)
    mle <- suppressWarnings(stats::glm.fit(mm, y, family = stats::binomial()))
    coef_mle <- mle$coefficients
    p_mle <- mle$fitted.values
    mle_stable <- isTRUE(mle$converged) && all(is.finite(coef_mle)) && max(abs(coef_mle)) < 15 && min(p_mle) >= 0.001 && 
        max(p_mle) <= 0.999
    if (mle_stable) {
        beta <- coef_mle
        pscore <- p_mle
        method <- "Unpenalized maximum-likelihood logistic regression"
        lambda <- 0
        convergence <- mle$converged
    }
    else {
        lambda <- ctx$ps_lambda
        ridge <- ridge_logistic(mm, y, lambda)
        beta <- ridge$par
        pscore <- plogis(drop(mm %*% beta))
        method <- paste0("Fixed-ridge logistic regression (stability fallback; lambda=", lambda, ")")
        convergence <- ridge$convergence == 0
    }
    pscore <- pmin(pmax(pscore, 1e-08), 1 - 1e-08)
    ow <- ifelse(y == 1, 1 - pscore, pscore)
    weights <- data.frame(Patient_Key = d$Patient_Key, Treatment_Group = y, Propensity_Score = pscore, Overlap_Weight = ow, 
        Analysis_Cohort = "Primary comparative cohort", Model_Method = method, stringsAsFactors = FALSE)
    coefficient <- data.frame(Term = colnames(mm), Coefficient = as_num(beta), Odds_Ratio_per_Model_Unit = exp(as_num(beta)), 
        Penalty_Lambda = lambda, stringsAsFactors = FALSE)
    centers <- data.frame(Covariate = c("Age", "Diagnosis_to_HSCT_Days_Calculated", "PreTx_Treatment_Cycles", "PreTx_BM_Blast_Pct", 
        "PreTx_FCM_Pct", "CD34_Dose_10e6_per_kg"), Imputation_Median = c(age$median, dx$median, cycles$median, bm$median, 
        fcm$median, cd34$median), Standardization_Center = c(z_age$center, z_dx$center, z_cycles$center, z_bm$center, z_fcm$center, 
        z_cd34$center), Standardization_Scale = c(z_age$scale, z_dx$scale, z_cycles$scale, z_bm$scale, z_fcm$scale, z_cd34$scale), 
        stringsAsFactors = FALSE)
    spec <- data.frame(Item = c("Estimand", "Cohort", "Outcome", "Model", "Formula", "Continuous missingness", "Categorical missingness", 
        "Overlap weights", "Variance in weighted outcome models", "Stability rule", "Unpenalized diagnostic", "Final penalty lambda", 
        "Reference-result policy", "Model convergence", "PS range"), Specification = c("Overlap-weighted average treatment effect in the overlap population", 
        "Primary comparative cohort: 47 treated and 23 platform-compatible controls", "Observed maintenance assignment", 
        method, formula_text, "Overall primary-cohort median imputation plus explicit missingness indicators", "Explicit Missing level if present", 
        "Treated: 1-PS; control: PS", "Robust sandwich variance clustered by Patient_Key", "Use unpenalized MLE only if converged, all coefficients finite, max |beta| <15, and fitted PS range is within 0.001-0.999; otherwise use fixed ridge lambda.", 
        paste0("converged=", mle$converged, "; max|beta|=", fmt_num(max(abs(coef_mle), na.rm = TRUE), 6), "; fitted PS range=", 
            fmt_num(min(p_mle, na.rm = TRUE), 8), " to ", fmt_num(max(p_mle, na.rm = TRUE), 8), "; stability_pass=", mle_stable), 
        as.character(lambda), "No tuning against reference results", as.character(convergence), paste0(fmt_num(min(pscore), 
            6), " to ", fmt_num(max(pscore), 6))), stringsAsFactors = FALSE)
    balance_rows <- list()
    treat <- y
    for (j in seq.int(2L, ncol(mm))) {
        x <- mm[, j]
        vt <- stats::var(x[treat == 1])
        vc <- stats::var(x[treat == 0])
        denom <- sqrt((vt + vc)/2)
        before <- smd_one(x, treat, rep(1, length(x)), denom)
        after <- smd_one(x, treat, ow, denom)
        balance_rows[[length(balance_rows) + 1L]] <- data.frame(Model_Term = colnames(mm)[j], Treated_Mean_Before = mean(x[treat == 
            1]), Control_Mean_Before = mean(x[treat == 0]), Treated_Mean_After = weighted_mean_safe(x[treat == 1], ow[treat == 
            1]), Control_Mean_After = weighted_mean_safe(x[treat == 0], ow[treat == 0]), SMD_Before = before, SMD_After = after, 
            Absolute_SMD_Before = abs(before), Absolute_SMD_After = abs(after), stringsAsFactors = FALSE)
    }
    balance <- do.call(rbind, balance_rows)
    ess_rows <- lapply(c(0, 1), function(g) {
        wg <- ow[y == g]
        data.frame(Group = if (g == 1) 
            "Treated"
        else "Control", N = length(wg), Sum_Weights = sum(wg), ESS = sum(wg)^2/sum(wg^2), Min_Weight = min(wg), Max_Weight = max(wg), 
            stringsAsFactors = FALSE)
    })
    ess <- do.call(rbind, ess_rows)
    ess <- rbind(ess, data.frame(Group = "Total", N = length(ow), Sum_Weights = sum(ow), ESS = sum(ow)^2/sum(ow^2), Min_Weight = min(ow), 
        Max_Weight = max(ow)))
    ps_dataset <- cbind(data.frame(Patient_Key = d$Patient_Key, Treatment_Group = y), psdf, weights[c("Propensity_Score", 
        "Overlap_Weight")])
    ctx$ps_weights <- weights
    ctx$ps_coefficients <- coefficient
    ctx$ps_centers <- centers
    ctx$ps_spec <- spec
    ctx$ps_balance <- balance
    ctx$ps_ess <- ess
    ctx$ps_dataset <- ps_dataset
    ctx$primary70$Propensity_Score <- pscore
    ctx$primary70$Overlap_Weight <- ow
    write_utf8_csv(spec, file.path(ctx$paths$output_data, "ps_model_specification.csv"))
    write_utf8_csv(coefficient, file.path(ctx$paths$output_data, "ps_model_coefficients.csv"))
    write_utf8_csv(centers, file.path(ctx$paths$output_data, "ps_imputation_standardization.csv"))
    write_utf8_csv(weights, file.path(ctx$paths$output_data, "patient_level_ps_overlap_weights.csv"))
    write_utf8_csv(ess, file.path(ctx$paths$output_data, "ess_overlap_weights.csv"))
    write_utf8_csv(balance, file.path(ctx$paths$output_data, "balance_smd_all_terms.csv"))
    write_utf8_csv(ps_dataset, file.path(ctx$paths$output_data, "ps_analysis_dataset.csv"))
    max_after <- max(balance$Absolute_SMD_After[is.finite(balance$Absolute_SMD_After)])
    add_validation(ctx, "RPL-PS-001", "Propensity score", "Single PS model", "One source-derived model; converged", paste(method, 
        "converged=", convergence), if (convergence) 
        "PASS"
    else "FAIL")
    add_validation(ctx, "RPL-PS-002", "Propensity score", "PS and overlap weights", "Finite PS in (0,1), positive finite weights", 
        paste0("PS ", fmt_num(min(pscore), 5), "-", fmt_num(max(pscore), 5), "; weight ", fmt_num(min(ow), 5), "-", fmt_num(max(ow), 
            5)), if (all(is.finite(pscore)) && all(pscore > 0 & pscore < 1) && all(is.finite(ow)) && all(ow > 0)) 
            "PASS"
        else "FAIL")
    add_validation(ctx, "RPL-PS-003", "Balance", "Model-matrix term maximum absolute SMD after weighting", "<=0.10", fmt_num(max_after, 
        6), if (max_after <= 0.1) 
        "PASS"
    else if (max_after <= 0.2) 
        "WARN"
    else "FAIL")
    add_validation(ctx, "RPL-PS-004", "Effective sample size", "Group ESS", "Both group ESS >5", paste(paste(ess$Group, fmt_num(ess$ESS, 
        3), sep = "="), collapse = "; "), if (all(ess$ESS[ess$Group != "Total"] > 5)) 
        "PASS"
    else "FAIL")
    append_log(ctx, paste0("Fit independent PS model: ", method, "; total ESS=", fmt_num(ess$ESS[ess$Group == "Total"], 3), 
        "; model-matrix term maximum |SMD| after=", fmt_num(max_after, 4), "."))
    ctx
}

pipeline_fit_survival_models <- 
function (ctx) 
{
    d <- ctx$primary70
    defs <- list(OS = c("OS_Months", "OS_Event"), LFS = c("LFS_Months", "LFS_Event"), GRFS_Main = c("GRFS_Months_Main", "GRFS_Event_Main"), 
        GRFS_Sensitivity = c("GRFS_Months_Sensitivity", "GRFS_Event_Sensitivity"))
    all_results <- list()
    all_intervals <- list()
    ctx$td_fits <- list()
    for (nm in names(defs)) {
        ans <- fit_td_endpoint(d, defs[[nm]][1], defs[[nm]][2], nm, ctx$days_per_month)
        all_results[[nm]] <- ans$result
        z <- ans$intervals
        z$Endpoint <- nm
        all_intervals[[nm]] <- z
        ctx$td_fits[[nm]] <- ans$fits
    }
    td_results <- do.call(rbind, all_results)
    td_intervals <- do.call(rbind, all_intervals)
    landmarks <- do.call(rbind, lapply(c(30, 60, 100), function(day) {
        do.call(rbind, lapply(c("OS", "LFS", "GRFS_Main"), function(nm) {
            fit_landmark(d, defs[[nm]][1], defs[[nm]][2], nm, day, ctx$days_per_month)
        }))
    }))
    ctx$td_results <- td_results
    ctx$td_intervals <- td_intervals
    ctx$landmark_results <- landmarks
    write_utf8_csv(td_results, file.path(ctx$paths$output_data, "time_dependent_cox_all_terms.csv"))
    write_utf8_csv(td_results[td_results$Term == "Exposure", , drop = FALSE], file.path(ctx$paths$output_data, "time_dependent_cox_exposure.csv"))
    write_utf8_csv(td_intervals, file.path(ctx$paths$output_data, "time_dependent_cox_intervals.csv"))
    write_utf8_csv(landmarks, file.path(ctx$paths$output_data, "landmark_sensitivity_results.csv"))
    core <- td_results$Term == "Exposure" & td_results$Model %in% c("Adjusted_unweighted", "Overlap_weighted", "Overlap_weighted_residual_adjusted")
    finite_core <- all(is.finite(td_results$HR[core])) && all(is.finite(td_results$CI_Lower_95[core])) && all(is.finite(td_results$CI_Upper_95[core]))
    add_validation(ctx, "RPL-TD-001", "Outcome model", "Primary time-dependent Cox models", "Exposure HR estimable for OS, LFS, main GRFS and sensitivity GRFS", 
        paste(sum(core), "exposure estimates; finite=", finite_core), if (sum(core) == 12L && finite_core) 
            "PASS"
        else "FAIL")
    add_validation(ctx, "RPL-TD-002", "Outcome model", "Delayed treatment initiation", "Exposure=0 until actual maintenance start and 1 thereafter", 
        paste0(nrow(td_intervals), " counting-process intervals"), if (all(td_intervals$tstop > td_intervals$tstart) && all(td_intervals$Exposure %in% 
            c(0, 1))) 
            "PASS"
        else "FAIL")
    add_validation(ctx, "RPL-TD-003", "Sensitivity", "Landmark analyses", "Day 30, 60 and 100 estimates for OS, LFS and main GRFS", 
        paste(nrow(landmarks), "term estimates"), if (nrow(landmarks) >= 9L) 
            "PASS"
        else "WARN")
    append_log(ctx, paste0("Fit adjusted and overlap-weighted time-dependent Cox models from HSCT origin, ", "plus prespecified landmark sensitivities. ", 
        sum(core), " core exposure estimates; finite=", finite_core, "."))
    ctx
}

ridge_logistic <- 
function (mm, y, lambda) 
{
    nll <- function(beta) {
        eta <- drop(mm %*% beta)
        p <- pmin(pmax(plogis(eta), 1e-12), 1 - 1e-12)
        -sum(y * log(p) + (1 - y) * log(1 - p)) + 0.5 * lambda * sum(beta[-1L]^2)
    }
    grad <- function(beta) {
        p <- plogis(drop(mm %*% beta))
        drop(crossprod(mm, p - y)) + c(0, lambda * beta[-1L])
    }
    stats::optim(rep(0, ncol(mm)), nll, grad, method = "BFGS", control = list(maxit = 2000, reltol = 1e-11))
}

smd_one <- 
function (x, treat, w = rep(1, length(x)), denominator = NULL) 
{
    x <- as_num(x)
    treat <- as.integer(treat)
    mt <- weighted_mean_safe(x[treat == 1], w[treat == 1])
    mc <- weighted_mean_safe(x[treat == 0], w[treat == 0])
    if (is.null(denominator)) {
        vt <- stats::var(x[treat == 1], na.rm = TRUE)
        vc <- stats::var(x[treat == 0], na.rm = TRUE)
        denominator <- sqrt((vt + vc)/2)
    }
    if (!is.finite(denominator) || denominator == 0) {
        return(if (isTRUE(all.equal(mt, mc))) 0 else Inf)
    }
    (mt - mc)/denominator
}

standardize_fixed <- 
function (x) 
{
    mu <- mean(x)
    s <- stats::sd(x)
    if (!is.finite(s) || s == 0) 
        s <- 1
    list(value = (x - mu)/s, center = mu, scale = s)
}

weighted_mean_safe <- 
function (x, w) 
{
    ok <- is.finite(x) & is.finite(w) & w >= 0
    if (!any(ok) || sum(w[ok]) <= 0) 
        return(NA_real_)
    sum(x[ok] * w[ok])/sum(w[ok])
}

weighted_mean_sd_text <- 
function (x, w) 
{
    m <- weighted_mean_safe(as_num(x), w)
    s <- sqrt(weighted_var_safe(as_num(x), w))
    paste0(fmt_num(m, 1), " (", fmt_num(s, 1), ")")
}

weighted_var_safe <- 
function (x, w) 
{
    ok <- is.finite(x) & is.finite(w) & w > 0
    if (sum(ok) < 2L) 
        return(NA_real_)
    mu <- weighted_mean_safe(x[ok], w[ok])
    sum(w[ok] * (x[ok] - mu)^2)/sum(w[ok])
}

write_utf8_csv <- 
function (x, path) 
{
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    write.csv(x, path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
    invisible(path)
}

