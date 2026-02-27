### fit_wrappers.R --- 
#----------------------------------------------------------------------
## Author: Johan Sebastian Ohlendorff
## Created: Feb 27 2026 (14:39) 
## Version: 
## Last-Updated: Feb 27 2026 (19:48) 
##           By: Johan Sebastian Ohlendorff
##     Update #: 31
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

hazard_fit <- function(data, model_hazard, outcome_string, covariates = NULL, formula_strategy = "additive", use_history_of_variables = FALSE, lag = NULL, k = NULL, time_covariates = NULL, baseline_covariates = NULL) {
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
            do.call(model_hazard, list(character_formula = formula_hazard, data = data))
        },
        error = function(e) {
            stop("Error in fitting hazard model: ", e, "with formula: ", formula_hazard)
        },
        warning = function(w) {
            message("Warning in fitting hazard model: ", w, "with formula: ", formula_hazard)
        }
        )
}

regression_fit <- function(data, model_regression, outcome_string, covariates = NULL, formula_strategy = "additive", use_history_of_variables = FALSE, lag = NULL, k = NULL, time_covariates = NULL, baseline_covariates = NULL, type = "propensity") {
       if (use_history_of_variables) {
           covariates <- get_history_of_variables(
                data,
                time_covariates,
                baseline_covariates,
                type = type,
                lag = lag,
                k = k
            )
       }
       
        withCallingHandlers(
        {
            if (type == "pseudo_outcome") {
                learn_Q(
                    model_type = model_regression,
                    history_of_variables = covariates,
                    data_learn = data,
                    formula_strategy = formula_strategy,
                    outcome_name = outcome_string
                )
            } 
        },
        error = function(e) {
            stop("Error in fitting regression model: ", e, "with formula: ", formula_regression)
        },
        warning = function(w) {
            message("Warning in fitting regression model: ", w, "with formula: ", formula_regression)
        }
        )
}

######################################################################
### fit_wrappers.R ends here
