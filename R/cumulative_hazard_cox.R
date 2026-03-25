### cumulative_hazard_cox.R --- 
#----------------------------------------------------------------------
## Author: Johan Sebastian Ohlendorff
## Created: Mar  4 2026 (16:29) 
## Version: 
## Last-Updated: Mar 25 2026 (11:47) 
##           By: Johan Sebastian Ohlendorff
##     Update #: 119
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

cumulative_hazard_cox <- function(fit, data, covariate_data, time_variable = "time",  time_ref = NULL){
    hazard <- hazard_minus <- Lambda <- exp_lp <- Lambda_minus <- hazard_prev <- exact_match <- hazard_time_ref <- NULL
    fit$coefficients[is.na(fit$coefficients)] <- 0 ## Force Brice's functions to behave
    ## Find exp(LP); i.e., exponential of linear predictor
    exp_lp_dt <- data.table(id = covariate_data$id)
    data.table::set(exp_lp_dt, j = "exp_lp", value = exp(riskRegression::coxLP(fit,data = covariate_data,center = FALSE)))

    ## Baseline cumulative hazard Lambda_0^x (T_j) for all j
    base_hazard <- riskRegression::predictCox(fit,centered = FALSE, type = "cumhazard")
    base_hazard <- data.table(hazard = base_hazard$cumhazard, time = base_hazard$time)
    setnames(base_hazard, "time", time_variable)
    N <- nrow(base_hazard)
    data.table::set(base_hazard, j = "hazard_prev", value = c(0, base_hazard$hazard[-N]))

    if (!is.null(time_ref)) {
        min_time <- min(data[, time_ref, with = FALSE])
        data_time_ref <- unique(data[, c("id", time_ref), with = FALSE])
        setnames(data_time_ref, time_ref, time_variable)
    } else {
        min_time <- 0
    }
    ## Add 0 row to base_hzazrd at min_time
    if (min_time < base_hazard[[time_variable]][1]) {
        temp_dat <- data.table(hazard = 0, time = min_time, hazard_prev = 0)
        setnames(temp_dat, "time", time_variable)
        base_hazard <- rbindlist(list(temp_dat, base_hazard))
    }
    matched_data <- base_hazard[data, on = time_variable, roll = TRUE]
    exact_matches <- which(data[[time_variable]] %in% base_hazard[[time_variable]])
    not_exact_matches <- setdiff(seq_len(nrow(matched_data)), exact_matches)
    data.table::set(matched_data, i = exact_matches, j = "hazard_minus", value = matched_data$hazard_prev[exact_matches])
    data.table::set(matched_data, i = not_exact_matches, j = "hazard_minus", value = matched_data$hazard[not_exact_matches])
    data.table::set(matched_data, j = "hazard_prev", value = NULL)
    if (!is.null(time_ref)){
        matched_data_time_ref <- base_hazard[data_time_ref, on = time_variable, roll = TRUE]
        ids_hazard_zero <- which(matched_data_time_ref[[time_variable]] == 0)
        set(matched_data_time_ref, i = ids_hazard_zero, j = "hazard", value = 0)
        setnames(matched_data_time_ref, "hazard", "hazard_time_ref")
        matched_data <- matched_data[matched_data_time_ref[, c("id", "hazard_time_ref"), with = FALSE], on = "id"]
    } else {
        set(matched_data, j = "hazard_time_ref", value = 0)
    }
        
    matched_data <- matched_data[
        exp_lp_dt,
        on = "id"
    ]
    set(
        matched_data,
        j = "Lambda",
        value = matched_data$exp_lp * (matched_data$hazard - matched_data$hazard_time_ref)
    )
    set(
        matched_data,
        j = "Lambda_minus",
        value = matched_data$exp_lp * (matched_data$hazard_minus - matched_data$hazard_time_ref)
    )
    set(
        matched_data,
        j = c("hazard", "hazard_minus", "hazard_time_ref"),
        value = NULL
    )
    matched_data
}

######################################################################
### cumulative_hazard_cox.R ends here
