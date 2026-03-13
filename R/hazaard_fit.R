### hazaard_fit.R --- 
#----------------------------------------------------------------------
## Author: Johan Sebastian Ohlendorff
## Created: Mar 13 2026 (18:42) 
## Version: 
## Last-Updated: Mar 13 2026 (18:42) 
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
hazard_fit <- function(data, model_hazard, outcome_string, covariates = NULL, formula_strategy = "additive", use_history_of_variables = FALSE, lag = NULL, k = NULL, time_covariates = NULL, baseline_covariates = NULL, time_variable = "time",penalize) {
       if (use_history_of_variables) {
           covariates <- get_history_of_variables(
                data,
                time_covariates,
                baseline_covariates,
                type = "hazard",
                lag = lag,
                k = k
            )
       }
       if (formula_strategy == "additive") {
            formula_hazard <- paste0(
                outcome_string, " ~ ",
                paste(covariates, collapse = "+")
            )
       } else {
           stop("Currently only 'additive' formula strategy is supported.")
       }
        withCallingHandlers(
        {
            do.call(model_hazard, list(character_formula = formula_hazard, data = data, time_variable = time_variable, penalize = penalize))
        },
        error = function(e) {
            stop("Error in fitting hazard model: ", e, "with formula: ", formula_hazard)
        },
        warning = function(w) {
            message("Warning in fitting hazard model: ", w, "with formula: ", formula_hazard)
        }
        )
}

######################################################################
### hazaard_fit.R ends here
