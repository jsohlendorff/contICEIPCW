### ipcw.R --- 
#----------------------------------------------------------------------
## Author: Johan Sebastian Ohlendorff
## Created: Mar  4 2026 (22:54) 
## Version: 
## Last-Updated: Mar  6 2026 (11:13) 
##           By: Johan Sebastian Ohlendorff
##     Update #: 43
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
        
        data_use[, c("time", "time_prev") := list(time_k, time_k_prev)]
        data_use[, type := "Gminus"]
        data_time_horizon <- copy(data_use)
        data_time_horizon[, type := "Gtau"]
        data_time_horizon[,time := time_horizon]
        dt <- rbind(data_use, data_time_horizon)

        data_use <- cumulative_hazard_cox(marginal_censoring_fit$fit, dt, data_use, time_ref = "time_prev")
        data_tau <- data_use[type == "Gtau"]
        data_tau[, Gtau := exp(-Lambda)]

        data_minus <- data_use[type == "Gminus"]
        data_minus[, Gminus := exp(-Lambda_minus)]
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
