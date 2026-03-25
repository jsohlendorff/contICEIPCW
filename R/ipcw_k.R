### ipcw_k.R --- 
#----------------------------------------------------------------------
## Author: Johan Sebastian Ohlendorff
## Created: Mar  4 2026 (22:54) 
## Version: 
## Last-Updated: Mar 25 2026 (14:44) 
##           By: Johan Sebastian Ohlendorff
##     Update #: 70
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

ipcw_k <- function(data, k, marginal_censoring_fit, time_horizon, is_censored, fast_ipcw = FALSE, survival_function) {
    ## FIXME: error when the pseudo outcome model is a weighted glm and marginal censoring is not assumed.
    if (is.null(marginal_censoring_fit)) {
        fast_ipcw <- TRUE
    }
    event_prev <- time_k <- time_k_prev <- NULL
    if (!is_censored) {
        Gminus <- Gtau <- rep(1, nrow(data))
    } else if (!fast_ipcw) {
        data_use <- data[event_k_prev %in% c("A", "L")]
        set(
            data_use,
            j = c("time", "time_prev"),
            value = list(
                data_use[[paste0("time_", k)]],
                data_use[[paste0("time_", k - 1)]]
            )
        )
        set(
            data_use,
            j = "type",
            value = "Gminus"
        )

        data_time_horizon <- copy(data_use)
        set(
            data_time_horizon,
            j = c("type","time"),
            value = list("Gtau", time_horizon)
        )

        dt <- rbind(data_use, data_time_horizon)

        data_use <- cumulative_hazard_cox(marginal_censoring_fit$fit,
                                          dt,
                                          data_use,
                                          time_ref = "time_prev",
                                          baseline_hazard = marginal_censoring_fit$baseline_hazard)
        data_tau <- data_use[type == "Gtau"]
        set(
            data_tau,
            j = "Gtau",
            value = exp(-data_tau$Lambda)
        )

        data_minus <- data_use[type == "Gminus"]
        set(
            data_minus,
            j = "Gminus",
            value = exp(-data_minus$Lambda_minus)
        )
        Gminus <- data_minus$Gminus
        Gtau <- data_tau$Gtau
    }
    event_k <- data[[paste0("event_", k)]]
    time_k <- data[[paste0("time_", k)]]
    if (!fast_ipcw) {
        ipcw <- 1*(event_k != "C" & time_k <= time_horizon) / Gminus + 1*(time_k > time_horizon) / Gtau
    } else {
        ipcw <- 1*(event_k != "C" & time_k <= time_horizon) / survival_function ## Only need Gtau for certain estimating equations anyway
    }
    
    return(ipcw)
}
######################################################################
### ipcw.R ends here
