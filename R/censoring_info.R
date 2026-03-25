### censoring_info.R --- 
#----------------------------------------------------------------------
## Author: Johan Sebastian Ohlendorff
## Created: Feb 27 2026 (14:15) 
## Version: 
## Last-Updated: Mar 25 2026 (13:07) 
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
censoring_info <- function(timevarying_data, baseline_data, time_horizons, marginal_censoring) {
    is_censored <- lapply(time_horizons, function(x) timevarying_data[event == "C" & time < x, .N] > 0)
    
    ## If marginal_censoring is TRUE, get data with time-varying covariates
    if (marginal_censoring && any(unlist(is_censored))) {
        data_marginal_censoring <-
            timevarying_data[event %in% c("tauend", "C", "Y", "D")][
                baseline_data, on = "id"
            ]
    } else {
        data_marginal_censoring <- NULL
    }
    list(is_censored = is_censored, data_marginal_censoring = data_marginal_censoring)
}

######################################################################
### censoring_info.R ends here

