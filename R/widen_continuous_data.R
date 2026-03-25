### widen_continuous_data.R --- 
#----------------------------------------------------------------------
## Author: Johan Sebastian Ohlendorff
## Created: Mar 13 2026 (18:51) 
## Version: 
## Last-Updated: Mar 25 2026 (12:55) 
##           By: Johan Sebastian Ohlendorff
##     Update #: 5
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
    set(
        data_wide,
        j = c("event_0", "time_0"),
        value = list("A", 0)
    )
    data_wide
}

######################################################################
### widen_continuous_data.R ends here
