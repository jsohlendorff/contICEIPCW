### cumulative_hazard_cox.R --- 
#----------------------------------------------------------------------
## Author: Johan Sebastian Ohlendorff
## Created: Mar  4 2026 (16:29) 
## Version: 
## Last-Updated: Mar  5 2026 (16:50) 
##           By: Johan Sebastian Ohlendorff
##     Update #: 77
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

cumulative_hazard_cox <- function(fit, data, covariate_data, time_variable = "time",  time_ref = NULL){
    hazard <- hazard_minus <- Lambda <- exp_lp <- Lambda_minus <- NULL
    fit$coefficients[is.na(fit$coefficients)] <- 0 ## Force Brice's functions to behave
    ## Find exp(LP); i.e., exponential of linear predictor
    exp_lp_dt <- data.table(id = covariate_data$id)
    exp_lp_dt$exp_lp <- exp(riskRegression::coxLP(fit,data = covariate_data,center = FALSE))

    ## Baseline cumulative hazard Lambda_0^x (T_j) for all j
    base_hazard <- riskRegression::predictCox(fit,centered = FALSE, type = "cumhazard")
    base_hazard <- data.table(hazard = base_hazard$cumhazard, time = base_hazard$time)
    setnames(base_hazard, "time", time_variable)
    base_hazard[, hazard_prev := c(0, hazard[-.N])]

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
    matched_data[, exact_match := data[[time_variable]] %in% base_hazard[[time_variable]]]
    matched_data[exact_match == FALSE, hazard_minus := hazard]
    matched_data[exact_match == TRUE, hazard_minus := hazard_prev]
    matched_data[, exact_match := NULL]
    matched_data[, hazard_prev := NULL]
    if (!is.null(time_ref)){
        matched_data_time_ref <- base_hazard[data_time_ref, on = time_variable, roll = TRUE]
        matched_data_time_ref[get(time_variable) == 0, "hazard" := 0]
        setnames(matched_data_time_ref, "hazard", "hazard_time_ref")
        matched_data <- merge(matched_data, matched_data_time_ref[, c("id", "hazard_time_ref"), with = FALSE],  by = "id", all.x = TRUE)
    } else {
        matched_data[, hazard_time_ref := 0]
    }
        
    matched_data <- merge(matched_data, exp_lp_dt, by = "id")
    matched_data[, Lambda := exp_lp * (hazard - hazard_time_ref)]
    matched_data[, Lambda_minus := exp_lp * (hazard_minus - hazard_time_ref)]
    matched_data[, c("exp_lp", "hazard", "hazard_minus", "hazard_time_ref") := NULL]
    matched_data
}

######################################################################
### cumulative_hazard_cox.R ends here
