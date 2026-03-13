### get_history_of_variables.R --- 
#----------------------------------------------------------------------
## Author: Johan Sebastian Ohlendorff
## Created: Mar 13 2026 (18:51) 
## Version: 
## Last-Updated: Mar 13 2026 (18:51) 
##           By: Johan Sebastian Ohlendorff
##     Update #: 1
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
                                     k) {
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
        ## Allow for A and L to occur at the same time
        time_history <- c(time_history, paste0("time_", k), paste0(setdiff(time_covariates, "A"), "_", k))
    } else if (type == "martingale") {
        time_history <- setdiff(time_history, paste0(setdiff(time_covariates, "A"), "_", k - 1))
    }

    ## Full history of variables, i.e., covariates used in regressions
    history_of_variables <- c(time_history, baseline_covariates)

    ## Remove variables from history_of_variables that do not have more than one value
    ## in the data
    if (!type == "martingale") {
        history_of_variables <- setdiff(
            history_of_variables,
            names(which(vapply(data[, .SD, .SDcols = history_of_variables], function(x) length(unique(x)) <= 1, FUN.VALUE = logical(1))))
        )
    }
    history_of_variables
}

######################################################################
### get_history_of_variables.R ends here
