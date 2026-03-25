### widen_continuous_data.R --- 
#----------------------------------------------------------------------
## Author: Johan Sebastian Ohlendorff
## Created: Mar 13 2026 (18:51) 
## Version: 
## Last-Updated: Mar 25 2026 (11:58) 
##           By: Johan Sebastian Ohlendorff
##     Update #: 3
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

# Function to widen continuous data from the long format to the wide format
widen_continuous_data <- function(timevarying_data, baseline_data, time_covariates) {
    data_wide <- data.table::dcast(timevarying_data,
                                   id ~ event_number,
                                   value.var = c("time", "event", time_covariates)
                                   )

    ## Merge with baseline data
    data_wide <- data_wide[baseline_data, on = "id"]
    data_wide[, c("event_0", "time_0") := list("A", 0)]
}

######################################################################
### widen_continuous_data.R ends here
