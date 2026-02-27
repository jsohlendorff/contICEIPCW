## Utility functions for the package

# Function to widen continuous data from the long format to the wide format
widen_continuous_data <- function(timevarying_data, baseline_data, time_covariates) {
    data_wide <- data.table::dcast(timevarying_data,
                                   id ~ event_number,
                                   value.var = c("time", "event", time_covariates)
                                   )

    ## Merge with baseline data
    data_wide <- merge(data_wide, baseline_data, by = "id")
    data_wide[, c("event_0", "time_0") := list("A", 0)]
}

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
            names(which(sapply(data[, ..history_of_variables], function(x) length(unique(x)) <= 1)))
        )
    }
    history_of_variables
}

get_at_risk_data <- function(data,
                             k,
                             time_horizon = NULL) {
    ## Remove unused factor levels
    if (k == 1) {
        at_risk_interevent <- at_risk_before_time_horizon <- data
    } else {
        at_risk_interevent <- data[event_k_previous %in% c("A", "L"), env = list(event_k_previous = paste0("event_", k - 1))]
        if (nrow(at_risk_interevent) == 0) {
            return(NULL)
        }
        at_risk_before_time_horizon <- at_risk_interevent[time_previous < time_horizon, env = list(time_previous = paste0("time_", k - 1))]
        if (nrow(at_risk_before_time_horizon) == 0) {
            at_risk_before_time_horizon <- NULL
        }
        
        ## Shift the interevent times according to time_(k-1); makes modeling more natural
        at_risk_interevent[, paste0("time_", k) := time_k - time_k_minus_1,
                           env = list(
                               time_k = paste0("time_", k),
                               time_k_minus_1 = paste0("time_", k - 1)
                           )
                           ]
        for (j in seq_len(k - 1)) {
            at_risk_interevent[, paste0("time_", j) := time_k_minus_1 - time_j,
                               env = list(
                                   time_k_minus_1 = paste0("time_", k - 1),
                                   time_j = paste0("time_", j)
                               )
                               ]
            at_risk_interevent[, paste0("event_", j) := droplevels(event_j),
                               env = list(
                                   event_j = paste0("event_", j)
                               )
                               ]
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
