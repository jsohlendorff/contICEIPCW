### select_event.R --- 
#----------------------------------------------------------------------
## Author: Johan Sebastian Ohlendorff
## Created: Feb 27 2026 (12:26) 
## Version: 
## Last-Updated: Feb 27 2026 (13:53) 
##           By: Johan Sebastian Ohlendorff
##     Update #: 14
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## Adaptively select last event based on the data if not provided
select_last_event <- function(timevarying_data, time_horizon, last_event) {
    if (is.null(last_event)) {
        at_risk_table <- timevarying_data[time < time_horizon & event %in% c("A", "L"), .N, by = "event_number"]
        if (nrow(at_risk_table) == 0) {
            max_event_number <- 1
            last_event <- 0
        } else {
            max_event_number <- max(at_risk_table$event_number)
            last_event <- at_risk_table[N > 40, event_number[.N]]
            if (last_event < max_event_number) {
                message(
                    "Adaptively selecting last event number (N <= 40). Event number: ",
                    last_event
                )
            }
        }
    }

    timevarying_data <- timevarying_data[event_number <= last_event | !(event %in% c("A", "L"))]
    timevarying_data <- timevarying_data[, event_number := seq_len(.N), by = id]
    return(last_event + 1)
}

######################################################################
### select_event.R ends here
