### predict_intervention.R --- 
#----------------------------------------------------------------------
## Author: Johan Sebastian Ohlendorff
## Created: Mar 13 2026 (18:50) 
## Version: 
## Last-Updated: Mar 13 2026 (18:50) 
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

## Wrapper function to predict the outcome under an intervention
predict_intervention <- function(data, k, predict_fun, static_intervention) {
  event_k <- A_0 <- event_k_intervention <- NULL
  intervened_data <- copy(data)
  if (k > 0) {
    intervened_data[event_k_intervention == "A", paste0("A_", k) := static_intervention, env = list(event_k_intervention = paste0("event_", k))]
  } else {
    intervened_data[, A_0 := static_intervention]
  }
  f <- predict_fun(intervened_data)

  ## Check if the predictions are in the range [0,1] if so warn and truncate
  if (any(f < 0 | f > 1)) {
    message("Predictions contain values outside the range [0, 1]. Truncating to [0, 1].")
    f <- pmin(pmax(f, 0), 1)
  }
  
  ## Warn if any predictions are NA or below or above 1
  if (any(is.na(f))) {
    stop("Predictions contain NA values.")
  }
  f
}

######################################################################
### predict_intervention.R ends here
