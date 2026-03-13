### learn_h2o.R --- 
#----------------------------------------------------------------------
## Author: Johan Sebastian Ohlendorff
## Created: Mar 13 2026 (18:49) 
## Version: 
## Last-Updated: Mar 13 2026 (19:06) 
##           By: Johan Sebastian Ohlendorff
##     Update #: 2
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
## Example function
## Export
##' Learn a model using H2O AutoML and return a prediction function
##' @param character_formula A character string representing the formula for the model, e.g., "outcome ~ var1 + var2".
##' @param data A data.table containing the data to learn the model from.
##' @param intervened_data A data.table containing the data to predict on under the intervention. If NULL, predictions will be made on the original data.
##' @param max_runtime_secs Maximum runtime for H2O AutoML in seconds.
##' @param nfolds Number of folds for cross-validation in H2O AutoML.
##' @param verbose Whether to print H2O AutoML output.
##' @param ... Additional arguments to pass to H2O AutoML.
##' @export 
learn_h2o <- function(character_formula,
                      data,
                      intervened_data = NULL,
                      max_runtime_secs = 30,
                      nfolds = 10,
                      verbose = FALSE,
                      ...) {
  requireNamespace("h2o", quietly = TRUE)
  formula_object <- as.formula(character_formula)
  outcome_name <- as.character(formula_object[[2]])
  history_of_variables <- labels(stats::terms(formula_object))
  data <- data[, c(outcome_name, history_of_variables), with = FALSE]
  ## Check if only 0/1 values in outcome_name
  if (all(data[[outcome_name]] %in% c(0, 1))) {
    distribution <- "bernoulli"
    data[[outcome_name]] <- as.factor(data[[outcome_name]]) # Convert to factor for classification
  } else {
    distribution <- "AUTO"
  }

  if (!verbose) sink("/dev/null") # Suppress H2O output
  suppressWarnings({
    h2o::h2o.init()
  })
  data_h2o <- h2o::as.h2o(data)

  ## AutoML
  aml <- h2o::h2o.automl(
    y = outcome_name,
    training_frame = data_h2o,
    max_runtime_secs = max_runtime_secs,
    nfolds = nfolds,
    distribution = distribution,
    sort_metric = "MSE",
    ...
  )
  if (!verbose) sink()

  best_model <- aml@leader
  leaderboard <- as.data.table(aml@leaderboard)

  ## Print leader and glm
  leader <- aml@leaderboard[1, ]
  lm_models <- leaderboard[grepl("LM", leaderboard$model_id), ]
  leaderboard <- rbind(leader, lm_models)
  print(leaderboard)

  if (distribution == "bernoulli") {
    ## For binary, we need to convert the predictions to a vector
    predict_fun <- function(data) {
      newdata_h2o <- h2o::as.h2o(data)
      as.vector(h2o::h2o.predict(best_model, newdata = newdata_h2o)$p1) # p1 for class 1 probability
    }
  } else {
    ## For regression, we can directly use the predict method
    predict_fun <- function(data) {
      newdata_h2o <- h2o::as.h2o(data)
      as.vector(h2o::h2o.predict(best_model, newdata = newdata_h2o)$predict)
    }
  }
  if (is.null(intervened_data)) {
    return(list(pred = predict_fun(data), predict_fun = predict_fun))
  } else {
    return(predict_fun(intervened_data))
  }
}

######################################################################
### learn_h2o.R ends here
