### get_history_of_variables.R --- 
#----------------------------------------------------------------------
## Author: Johan Sebastian Ohlendorff
## Created: Mar 13 2026 (18:51) 
## Version: 
## Last-Updated: Mar 19 2026 (10:12) 
##           By: Johan Sebastian Ohlendorff
##     Update #: 15
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

get_history_of_variables <- function(data,
                                     time_covariates,
                                     baseline_covariates,
                                     type,
                                     lag,
                                     k,
                                     exclude_latest_covariate) {
    if (!is.null(lag) && k > 1) {
        event_points <- seq(from = max(1, k - lag), to = k - 1, by = 1)
    } else {
        event_points <- seq_len(k - 1)
    }

    ## Time-varying covariates to use in regressions
    if (k > 1) {
        time_history <- unlist(lapply(c(time_covariates, "time", "event"), function(x) {
            paste0(x, "_", event_points)
        }))
    } else {
        time_history <- NULL
    }
    if (type == "hazard") {
        time_history <- setdiff(time_history, paste0("time_", k - 1))
    } else if (type == "propensity") {
        time_covariates <- c("time", time_covariates)
        if (is.null(exclude_latest_covariate)) {
            time_history <- c(time_history, paste0(setdiff(time_covariates, "A"), "_", k))
        } else {
            time_history <- c(time_history, paste0(setdiff(time_covariates, c("A", exclude_latest_covariate)), "_", k))
        }
    } else if (type == "martingale") {
        time_history <- setdiff(time_history, paste0(setdiff(time_covariates, "A"), "_", k - 1))
    }

    ## Full history of variables, i.e., covariates used in regressions
    history_of_variables <- c(time_history, baseline_covariates)
    setdiff(
            history_of_variables,
        names(which(vapply(data[, .SD, .SDcols = history_of_variables], function(x) length(unique(x)) <= 1, FUN.VALUE = logical(1))))
    )
}

######################################################################
### get_history_of_variables.R ends here
