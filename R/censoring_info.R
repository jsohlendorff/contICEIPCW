### censoring_info.R --- 
#----------------------------------------------------------------------
## Author: Johan Sebastian Ohlendorff
## Created: Feb 27 2026 (14:15) 
## Version: 
## Last-Updated: Feb 27 2026 (16:33) 
##           By: Johan Sebastian Ohlendorff
##     Update #: 10
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
censoring_info <- function(timevarying_data, baseline_data, time_horizon, marginal_censoring, model_hazard, model_pseudo_outcome) {
    is_censored <- timevarying_data[event == "C" & time < time_horizon, .N] > 0

    ## If marginal_censoring is TRUE, get data with time-varying covariates
    if (marginal_censoring) {
        data_baseline <- merge(timevarying_data[event %in% c("tauend", "C", "Y", "D")],
                                baseline_data,
                                by = "id",
                                all.x = TRUE)
    } else {
        data_baseline <- NULL
    }
    
    ## Check user input if censored
    if (is_censored && is.null(model_hazard)) {
        stop(
            "Censoring is present, but no censoring model is provided.
             Please provide a censoring model such as `model_hazard = 'learn_coxph'`."
        )
    }

    if (model_pseudo_outcome == "quasibinomial" && is_censored) {
        stop("quasibinomial model is not suitable for censored data.
              Use scaled_quasibinomial instead.")
    }
    list(is_censored = is_censored, data_baseline = data_baseline)
}

######################################################################
### censoring_info.R ends here
