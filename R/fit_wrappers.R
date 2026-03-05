### fit_wrappers.R --- 
#----------------------------------------------------------------------
## Author: Johan Sebastian Ohlendorff
## Created: Feb 27 2026 (14:39) 
## Version: 
## Last-Updated: Mar  5 2026 (15:10) 
##           By: Johan Sebastian Ohlendorff
##     Update #: 53
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

regression_fit <- function(data, model_regression, outcome_string, outcome_string_unweighted = NULL, ipcw_name = NULL, covariates = NULL, formula_strategy = "additive", use_history_of_variables = FALSE, lag = NULL, k = NULL, time_covariates = NULL, baseline_covariates = NULL, type = "propensity", penalize) {
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
                    outcome_name = outcome_string,
                    outcome_string_unweighted = outcome_string_unweighted,
                    ipcw_name = ipcw_name,
                    penalize
                )
            } else if (type == "propensity") {
                formula_propensity <- paste0(
                    outcome_string, " ~ ",
                    paste(covariates, collapse = "+")
                )
                do.call(model_regression, list(character_formula = formula_propensity, data = data, penalize = penalize))
            } else {
                stop("Unsupported regression type: ", type)
            }
        },
        error = function(e) {
            stop("Error in fitting regression model: ", e, "with outcome: ", outcome_string, " and covariates: ", paste(covariates, collapse = ", "), " and type: ", type)
        },
        warning = function(w) {
            message("Warning in fitting regression model: ", w, "with outcome: ", outcome_string, " and covariates: ", paste(covariates, collapse = ", "), " and type: ", type)
        }
        )
}

######################################################################
### fit_wrappers.R ends here
