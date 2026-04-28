### predict_intervention.R --- 
#----------------------------------------------------------------------
## Author: Johan Sebastian Ohlendorff
## Created: Mar 13 2026 (18:50) 
## Version: 
## Last-Updated: Apr 28 2026 (12:08) 
##           By: Johan Sebastian Ohlendorff
##     Update #: 26
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## Wrapper function to predict the outcome under an intervention
predict_intervention <- function(data, k, predict_fun, static_intervention, verbose, intervention_hook, intervene_all_sequential_regression) {
  intervened_data <- copy(data)
  if (!intervene_all_sequential_regression) {
      if (k > 0) {
          which_event_A <- which(intervened_data[[paste0("event_", k)]] == "A")
          set(
              intervened_data,
              i = which_event_A,
              j = paste0("A_", k),
              value = static_intervention
          )
      } else {
          set(
              intervened_data,
              j = "A_0",
              value = static_intervention
          )
      } 
  } else {
      for (j in 0:k){
          set(
              intervened_data,
              j = paste0("A_", j),
              value = static_intervention
          )
      }
  }
  
  if (!is.null(intervention_hook)) {
    intervened_data <- intervention_hook(intervened_data, k)
  }
  f <- predict_fun(intervened_data)

  ## Check if the predictions are in the range [0,1] if so warn and truncate
  if (any(f < 0 | f > 1)) {
    if (verbose) {
      message("Predictions contain values outside the range [0, 1]. Truncating to [0, 1].")
    }
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
