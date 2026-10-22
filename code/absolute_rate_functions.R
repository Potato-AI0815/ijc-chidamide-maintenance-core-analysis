# Existing current-cohort KM and CIF functions; no raw-data paths.

as_int <- 
function (x) 
suppressWarnings(as.integer(as.character(x)))

step_value <- 
function (time, value, horizon) 
{
    ok <- is.finite(time) & !is.na(value) & time <= horizon
    if (!any(ok)) 
        return(0)
    idx <- which(ok)
    max_time <- max(time[idx], na.rm = TRUE)
    value[max(idx[time[idx] == max_time], na.rm = TRUE)]
}

make_competing_status <- 
function (event_time, event_flag, competing_times, censor_time) 
{
    n <- length(event_flag)
    out_time <- as_num(censor_time)
    out_status <- rep(0L, n)
    event_time <- as_num(event_time)
    event_flag <- as_int(event_flag)
    competing_times <- as.matrix(competing_times)
    for (i in seq_len(n)) {
        target <- if (!is.na(event_flag[i]) && event_flag[i] == 1L && is.finite(event_time[i])) 
            event_time[i]
        else Inf
        comp <- as_num(competing_times[i, ])
        comp <- comp[is.finite(comp)]
        ctime <- if (length(comp)) 
            min(comp)
        else Inf
        censor <- if (is.finite(out_time[i])) 
            out_time[i]
        else Inf
        first <- min(target, ctime, censor)
        out_time[i] <- first
        if (is.finite(target) && target <= ctime && target <= censor) {
            out_status[i] <- 1L
        }
        else if (is.finite(ctime) && ctime < target && ctime <= censor) {
            out_status[i] <- 2L
        }
        else {
            out_status[i] <- 0L
        }
    }
    list(time = out_time, status = out_status)
}

km_points <- 
function (df, time_col, event_col, horizons, endpoint, group_label, notes) 
{
    time <- as_num(df[[time_col]])
    event <- as_int(df[[event_col]])
    ok <- is.finite(time) & !is.na(event)
    if (!any(ok)) 
        stop("No valid observations for KM endpoint: ", endpoint, " / ", group_label)
    fit <- survival::survfit(survival::Surv(time[ok], event[ok]) ~ 1, conf.type = "log-log")
    s <- summary(fit, times = horizons, extend = TRUE)
    out <- data.frame(endpoint = endpoint, time_point = paste0(horizons, " months"), group = group_label, N = sum(ok), events_or_target_events = vapply(horizons, 
        function(h) sum(event[ok] == 1L & time[ok] <= h, na.rm = TRUE), integer(1)), competing_events = NA_integer_, number_at_risk = s$n.risk, 
        censoring_count = vapply(horizons, function(h) sum(event[ok] == 0L & time[ok] <= h, na.rm = TRUE), integer(1)), estimate = s$surv, 
        lower_95_CI = s$lower, upper_95_CI = s$upper, method = "KM", notes = notes, stringsAsFactors = FALSE)
    out
}

cif_points <- 
function (time, status, horizons, endpoint, group_label, horizon_unit, notes, N_override = NULL) 
{
    time <- as_num(time)
    status <- as_int(status)
    ok <- is.finite(time) & !is.na(status)
    if (!any(ok)) 
        stop("No valid observations for CIF endpoint: ", endpoint, " / ", group_label)
    target_total <- sum(status[ok] == 1L)
    ci <- cmprsk::cuminc(ftime = time[ok], fstatus = status[ok], cencode = 0)
    curves <- ci[vapply(ci, is.list, logical(1))]
    curve <- NULL
    if ("1 1" %in% names(curves)) 
        curve <- curves[["1 1"]]
    if (is.null(curve)) {
        target_name <- grep("(^| )1$", names(curves), value = TRUE)
        if (length(target_name)) 
            curve <- curves[[target_name[1]]]
    }
    if (is.null(curve) && target_total > 0) 
        curve <- curves[[1L]]
    if (is.null(curve)) {
        curve <- list(time = numeric(0), est = numeric(0), var = numeric(0))
    }
    do.call(rbind, lapply(horizons, function(h) {
        est <- step_value(curve$time, curve$est, h)
        variance <- step_value(curve$time, curve$var, h)
        se <- sqrt(max(variance, 0))
        data.frame(endpoint = endpoint, time_point = if (horizon_unit == "days") 
            paste0(h, " days")
        else paste0(h, " months"), group = group_label, N = if (is.null(N_override)) 
            sum(ok)
        else N_override, events_or_target_events = sum(status[ok] == 1L & time[ok] <= h, na.rm = TRUE), competing_events = sum(status[ok] == 
            2L & time[ok] <= h, na.rm = TRUE), number_at_risk = NA_integer_, censoring_count = sum(status[ok] == 0L & time[ok] <= 
            h, na.rm = TRUE), estimate = est, lower_95_CI = max(0, est - 1.96 * se), upper_95_CI = min(1, est + 1.96 * se), 
            method = "CIF", notes = notes, stringsAsFactors = FALSE)
    }))
}

