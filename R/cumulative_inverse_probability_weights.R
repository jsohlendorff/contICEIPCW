### cumulative_inverse_probability_weights.R --- 
#----------------------------------------------------------------------
## Author: Johan Sebastian Ohlendorff
## Created: Feb 27 2026 (18:43) 
## Version: 
## Last-Updated: Mar 13 2026 (18:43) 
##           By: Johan Sebastian Ohlendorff
##     Update #: 28
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
## Estimate IPW weights in efficient influence function
cumulative_inverse_probability_weights <- function(data, static_intervention, time_horizon, return_ipw, last_event) {
    survival_censoring_0 <- cum_treatment_k <- cum_propensity_k <- cum_survival_censoring_k <- event_k_prev <- ipw_k <- event_k <- time_k <- survival_censoring_k <- ipw_cum_weight <- NULL
    ## Cumulative product of propensity scores (treatment)
    data[, paste0("propensity_", last_event) := 1]
    propensity_cols <- paste0("propensity_", seq(0, last_event))
    data[, (propensity_cols) := lapply(.SD, function(x) data.table::fifelse(is.na(x), 1, x)), .SDcols = propensity_cols]
    data[, paste0("cum_propensity_", seq(0, last_event)) := Reduce(`*`, .SD, accumulate = TRUE), .SDcols = propensity_cols]

    ## Cumulative product of censoring probabilities
    survival_censoring_cols <- paste0("survival_censoring_", seq(0, last_event))
    data[, survival_censoring_0 := 1]
    data[, paste0("cum_survival_censoring_", seq(0, last_event)) := Reduce(`*`, .SD, accumulate = TRUE), .SDcols = survival_censoring_cols]

    ## Cumulative product of whether treated according to static intervention is followed
    treatment_cols <- paste0("A_", seq(0, last_event - 1))
    data[, paste0("I_", seq(0, last_event - 1)) := lapply(.SD, function(col) as.integer(col == static_intervention)), .SDcols = treatment_cols]
    data[, paste0("cum_treatment_", seq(0, last_event - 1)) := Reduce(`*`, .SD, accumulate = TRUE), .SDcols = treatment_cols]
    data[, paste0("cum_treatment_", last_event) := 1]

    ## Calculate the inverse probability weights for the efficient influence function
    for (k in seq(0, last_event)) {
        data[, paste0("ipw_cum_weight_", k) := cum_treatment_k / (cum_propensity_k * cum_survival_censoring_k), env = list(                                                                      
            cum_treatment_k = paste0("cum_treatment_", k),
            cum_propensity_k = paste0("cum_propensity_", k),
            cum_survival_censoring_k = paste0("cum_survival_censoring_", k)
        )]
    }

    ## Calculate the inverse probability weights for the relevant estimator
    if (return_ipw) {
        for (k in seq(1, last_event)) {
            data[, paste0("ipw_", k) := 0]
            data[event_k_prev %in% c("A", "L"), ipw_k := (1 * (event_k == "Y" & time_k <= time_horizon)) / (survival_censoring_k) * ipw_cum_weight, env = list(
                                                                                                                            survival_censoring_k = paste0("survival_censoring_", k),
                                                                                                                            ipw_k = paste0("ipw_", k),
                                                                                                                            event_k = paste0("event_", k),
                                                                                                                            time_k = paste0("time_", k),
                                                                                                                            event_k_prev = paste0("event_", k - 1),
                                                                                                                            ipw_cum_weight = paste0("ipw_cum_weight_", k-1)
                                                                                                                        )]
        }
    }
}


######################################################################
### inverse_probability_weights.R ends here
