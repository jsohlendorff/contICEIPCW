### select_last_event.R --- 
#----------------------------------------------------------------------
## Author: Johan Sebastian Ohlendorff
## Created: Feb 27 2026 (12:26) 
## Version: 
## Last-Updated: Mar 31 2026 (11:30) 
##           By: Johan Sebastian Ohlendorff
##     Update #: 107
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## Adaptively select last event based on the data if not provided
select_last_event <- function(timevarying_data, time_horizons, last_non_terminal_event, min_events, verbose) {
    max_time_horizon <- max(time_horizons)
    if (is.null(last_non_terminal_event)) {
        last_non_terminal_event <- list()
        for (i in seq_len(length(time_horizons))) {
            at_risk_table <- timevarying_data[time < time_horizons[i] & event %in% c("A", "L"), .N, by = "event_number"]
            if (nrow(at_risk_table) == 0) {
                max_event_number <- 1
                last_non_terminal_event[[i]] <- 0
            } else {
                max_event_number <- max(at_risk_table$event_number)
                last_non_terminal_event[[i]] <- at_risk_table[N > min_events, event_number[.N]]
                if (last_non_terminal_event[[i]] < max_event_number && verbose) {
                    message(
                        "Adaptively selecting last event number (N <= ", min_events, "). Event number: ",
                        last_non_terminal_event[[i]]
                    )
                }
            }
        }
    }

    ## We should only start at the last event at which there is a terminal events, otherwise the iterative regression would just be regressing zero outcomes
    max_event_number_outcome <- lapply(time_horizons, function(x) timevarying_data[time <= x & event %in% "Y", .(max(event_number))]$V1)

    timevarying_data <- timevarying_data[event_number <= last_non_terminal_event[[which(max_time_horizon == time_horizons)]] | !(event %in% c("A", "L"))]
    timevarying_data <- timevarying_data[, event_number := seq_len(.N), by = id]
    
    info <- data.table::data.table(
        time_horizon = time_horizons,
        last_non_terminal_event = unlist(last_non_terminal_event),
        max_event_number_outcome = unlist(max_event_number_outcome)
    )
    info[, last_event := pmin(last_non_terminal_event + 1, max_event_number_outcome)]
    return(list(timevarying_data = timevarying_data,
                info = info))
}

######################################################################
### select_event.R ends here
