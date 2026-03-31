### regression_fit.R --- 
#----------------------------------------------------------------------
## Author: Johan Sebastian Ohlendorff
## Created: Mar 13 2026 (18:42) 
## Version: 
## Last-Updated: Mar 31 2026 (11:00) 
##           By: Johan Sebastian Ohlendorff
##     Update #: 84
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
regression_fit <- function(data,
                           model_regression,
                           outcome_string,
                           outcome_string_unweighted = NULL,
                           ipcw_name = NULL,
                           covariates = NULL,
                           formula_strategy = "additive",
                           use_history_of_variables = FALSE,
                           lag = NULL,
                           k = NULL,
                           time_covariates = NULL,
                           baseline_covariates = NULL,
                           type = "propensity",
                           penalize,
                           exclude_latest_covariate,
                           verbose,
                           reduce_colinearity_time = FALSE,
                           time_horizon = NULL) {
       if (use_history_of_variables) {
           covariates <- get_history_of_variables(
                data = data,
                time_covariates,
                baseline_covariates,
                type = type,
                lag = lag,
                k = k,
                exclude_latest_covariate = exclude_latest_covariate
            )
       }
       if (reduce_colinearity_time) {
           ## Avoid changing the original data by reference when modifying it for decolinearization
           ## FIXME: This is a bit hacky and may be slow.
           data <- copy(data)
           decolinearize_time <- function(data, k, type, time_horizon) {
               res <- list()
               if (type == "pseudo_outcome") {
                   if (k>1) {
                   res[[k-1]] <- time_horizon - data[[paste0("time_", k-1)]]
                   for (j in seq_len(k-2)) {
                       res[[j]] <- data[[paste0("time_", j)]] - data[[paste0("time_", j-1)]]
                   }
                   for (j in seq_len(k-1)) {
                       set(data, j = paste0("time_", j), value = res[[j]])
                   }
                   }
               } else {
                   for (j in seq_len(k)) {
                       res[[j]] <- data[[paste0("time_", j)]] - data[[paste0("time_", j-1)]]
                   }
                   for (j in seq_len(k)) {
                       set(data, j = paste0("time_", j), value = res[[j]])
                   }
               }
               return(data)
           }
           data <- decolinearize_time(data, k, type, time_horizon)
       }
       
        withCallingHandlers(
        {
            if (type == "pseudo_outcome") {
                fit <- learn_Q(
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
                fit <- do.call(model_regression, list(character_formula = formula_propensity, data = data, penalize = penalize))
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

       if (reduce_colinearity_time && type == "pseudo_outcome") {
           predict_fun <- function(newdata) {
               newdata <- copy(newdata)
               newdata <- decolinearize_time(newdata, k, type, time_horizon)
               fit(newdata)
             }
            return(predict_fun)
       } else {
            return(fit)
        }
}

######################################################################
### regression_fit.R ends here
