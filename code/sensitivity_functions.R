# Existing leave-one-comparator-out and alternative-analysis estimators; original weights retained.

extract_exposure <- 
function (fit, endpoint, model, n_patients, n_intervals, events, notes) 
{
    s <- summary(fit)
    co <- as.data.frame(s$coefficients)
    if (!("Exposure" %in% rownames(co))) 
        stop("Exposure term not found")
    se_col <- if ("robust se" %in% names(co)) 
        "robust se"
    else "se(coef)"
    p_col <- grep("^Pr\\(", names(co), value = TRUE)[1]
    beta <- co["Exposure", "coef"]
    se <- co["Exposure", se_col]
    data.frame(Endpoint = endpoint, Analysis = model, Term = "Exposure", Log_HR = beta, Robust_SE = se, HR = exp(beta), CI_Lower_95 = exp(beta - 
        1.96 * se), CI_Upper_95 = exp(beta + 1.96 * se), P_Value = co["Exposure", p_col], N_Patients = n_patients, N_Intervals = n_intervals, 
        Events = events, Notes = notes, stringsAsFactors = FALSE)
}

safe_fit <- 
function (expr) 
{
    warn <- character()
    val <- tryCatch(withCallingHandlers(expr, warning = function(w) {
        warn <<- c(warn, conditionMessage(w))
        invokeRestart("muffleWarning")
    }), error = function(e) e)
    list(value = val, warning = paste(unique(warn), collapse = " | "))
}

fit_weighted_td <- 
function (cp, endpoint, analysis, notes) 
{
    pats <- unique(cp$Patient_Key)
    fit <- survival::coxph(survival::Surv(tstart, tstop, event) ~ Exposure + survival::cluster(Patient_Key), data = cp, weights = Overlap_Weight, 
        ties = "breslow", robust = TRUE, model = TRUE, x = TRUE)
    extract_exposure(fit = fit, endpoint = endpoint, model = analysis, n_patients = length(pats), n_intervals = nrow(cp), 
        events = sum(cp$event == 1, na.rm = TRUE), notes = notes)
}

