### regression_fit.R --- 
#----------------------------------------------------------------------
## Author: Johan Sebastian Ohlendorff
## Created: Mar 13 2026 (18:42) 
## Version: 
## Last-Updated: Mar 18 2026 (14:50) 
##           By: Johan Sebastian Ohlendorff
##     Update #: 8
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
regression_fit <- function(data, model_regression, outcome_string, outcome_string_unweighted = NULL, ipcw_name = NULL, covariates = NULL, formula_strategy = "additive", use_history_of_variables = FALSE, lag = NULL, k = NULL, time_covariates = NULL, baseline_covariates = NULL, type = "propensity", penalize, verbose) {
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
                    penalize = penalize,
                    verbose = verbose,
                    k
                )
            } else if (type == "propensity") {
                formula_propensity <- paste0(
                    outcome_string, " ~ ",
                    paste(covariates, collapse = "+")
                )
                if (verbose) message("Fitting propensity score model with formula: ", formula_propensity)
                do.call(model_regression, list(character_formula = formula_propensity, data = data, penalize = penalize))
            } else {
                stop("Unsupported regression type: ", type)
            }
        },
        error = function(e) {
            stop("Error in fitting regression model: ", e, "with outcome: ", outcome_string, " and covariates: ", paste(covariates, collapse = ", "), " and type: ", type)
        },
        warning = function(w) {
            if (verbose) message("Warning in fitting regression model: ", w, "with outcome: ", outcome_string, " and covariates: ", paste(covariates, collapse = ", "), " and type: ", type)
        }
        )
}

######################################################################
### regression_fit.R ends here
