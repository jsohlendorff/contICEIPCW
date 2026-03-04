### censoring_info.R --- 
#----------------------------------------------------------------------
## Author: Johan Sebastian Ohlendorff
## Created: Feb 27 2026 (14:15) 
## Version: 
## Last-Updated: Mar  4 2026 (19:42) 
##           By: Johan Sebastian Ohlendorff
##     Update #: 17
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
censoring_info <- function(timevarying_data, baseline_data, time_horizon, marginal_censoring) {
    event<-time<-NULL
    is_censored <- timevarying_data[event == "C" & time < time_horizon, .N] > 0

    ## If marginal_censoring is TRUE, get data with time-varying covariates
    if (marginal_censoring) {
        data_marginal_censoring <- merge(timevarying_data[event %in% c("tauend", "C", "Y", "D")],
                                baseline_data,
                                by = "id",
                                all.x = TRUE)
    } else {
        data_marginal_censoring <- NULL
    }
    list(is_censored = is_censored, data_marginal_censoring = data_marginal_censoring)
}

######################################################################
### censoring_info.R ends here

