### cumulative_inverse_probability_weights.R --- 
#----------------------------------------------------------------------
## Author: Johan Sebastian Ohlendorff
## Created: Feb 27 2026 (18:43) 
## Version: 
## Last-Updated: Apr 28 2026 (11:57) 
##           By: Johan Sebastian Ohlendorff
##     Update #: 104
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
    ## Cumulative product of propensity scores (treatment)
    propensity_cols <- paste0("propensity_", seq(0, last_event - 1))
    data[, (propensity_cols) := lapply(.SD, function(x) data.table::fifelse(is.na(x), 1, x)), .SDcols = propensity_cols]
    data[, paste0("cum_propensity_", seq(0, last_event - 1)) := Reduce(`*`, .SD, accumulate = TRUE), .SDcols = propensity_cols]

    ## Cumulative product of censoring probabilities
    survival_censoring_cols <- paste0("survival_censoring_", seq(0, last_event))
    set(
        data,
        j = "survival_censoring_0",
        value = 1
    )
    data[, paste0("cum_survival_censoring_", seq(0, last_event)) := Reduce(`*`, .SD, accumulate = TRUE), .SDcols = survival_censoring_cols]

    ## Cumulative product of whether treated according to static intervention is followed
    treatment_cols <- paste0("A_", seq(0, last_event - 1))
    data[, paste0("I_", seq(0, last_event - 1)) := lapply(.SD, function(col) as.integer(col == static_intervention)), .SDcols = treatment_cols]
    data[, paste0("cum_treatment_", seq(0, last_event - 1)) := Reduce(`*`, .SD, accumulate = TRUE), .SDcols = treatment_cols]

    ## Calculate the inverse probability weights for the efficient influence function
    for (k in seq(0, last_event - 1)) {
        which_following <- which(data[[paste0("cum_treatment_", k)]] != 0)
        set(data,
            j = paste0("ipw_cum_weight_", k),
            value = 0
        )
        set(
            data,
            i = which_following,
            j = paste0("ipw_cum_weight_", k),
            value = data[[paste0("cum_treatment_", k)]][which_following] / (data[[paste0("cum_propensity_", k)]][which_following] * data[[paste0("cum_survival_censoring_", k)]][which_following])
        )
        ## Check if any weights are infinite
        if (any(is.infinite(data[[paste0("ipw_cum_weight_", k)]]))) {
            stop(paste0("Infinite weights detected at k=", k, " indicating zero probability of treatment among subjects following the regimen. You may want to consider truncating the weights by applying the argument gbounds in propensity_scores."))
        }
    }

    ## Calculate the inverse probability weights for the relevant estimator
    if (return_ipw) {
        for (k in seq(1, last_event)) {
            val <- (1 * (data[[paste0("event_", k)]] == "Y" & data[[paste0("time_", k)]] <= time_horizon)) / (data[[paste0("survival_censoring_", k)]]) * data[[paste0("ipw_cum_weight_", k-1)]]
            val[is.na(val)] <- 0
            set(data, j = paste0("ipw_", k), value = val)
        }
    }
}


######################################################################
### inverse_probability_weights.R ends here
