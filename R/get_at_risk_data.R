### get_at_risk_data.R --- 
#----------------------------------------------------------------------
## Author: Johan Sebastian Ohlendorff
## Created: Mar 13 2026 (18:52) 
## Version: 
## Last-Updated: Mar 13 2026 (18:52) 
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

get_at_risk_data <- function(data,
                             k,
                             time_horizon = NULL) {
    event_k_previous <- time_previous <- time_k <- time_k_prev_1 <- time_j <- event_j <- event_k <- time_k_prev <- event_k_prev <- NULL
    
    ## Create shortcuts for the k'th iteration
    data[, c("event_k", "time_k", "time_k_prev", "event_k_prev")
         := list(event_k, time_k, time_k_prev, event_k_prev),
         env = list(event_k = paste0("event_", k),
                    time_k = paste0("time_", k),
                    event_k_prev = paste0("event_", k - 1),
                    time_k_prev = paste0("time_", k - 1))]

    ## Remove unused factor levels
    if (k == 1) {
        at_risk_interevent <- at_risk_before_time_horizon <- data
    } else {
        at_risk_interevent <- data[event_k_prev %in% c("A", "L")]
        if (nrow(at_risk_interevent) == 0) {
            return(NULL)
        }
        at_risk_before_time_horizon <- at_risk_interevent[time_k_prev < time_horizon]
        if (nrow(at_risk_before_time_horizon) == 0) {
            at_risk_before_time_horizon <- NULL
        }
        
        ## Shift the interevent times according to time_(k-1); makes modeling more natural
        at_risk_interevent[, paste0("time_", k) := time_k - time_k_prev]
        for (j in seq_len(k - 1)) {
            at_risk_interevent[, paste0("time_", j) := time_k_prev - time_j,
                               env = list(time_j = paste0("time_", j))]
            at_risk_interevent[, paste0("event_", j) := droplevels(event_j),
                               env = list(event_j = paste0("event_", j))]
        }
    }
    at_risk_interevent[, (names(.SD)) := lapply(.SD, droplevels), 
                       .SDcols = is.factor]
    if (!is.null(at_risk_before_time_horizon)) {
        at_risk_before_time_horizon[, (names(.SD)) := lapply(.SD, droplevels), 
                                    .SDcols = is.factor]
    }
    list(
        at_risk_interevent = at_risk_interevent,
        at_risk_before_time_horizon = at_risk_before_time_horizon
    )
}

######################################################################
### get_at_risk_data.R ends here
