expit <- function(x) 1 / (1 + exp(-x))
logit <- function(p) log(p / (1 - p))
probit <- function(x) stats::pnorm(x)
inv_probit <- function(p) stats::qnorm(p)

# Model to use for the outcome regression which returns a prediction function
# for the chosen model.
# Available models are:
# \code{"quasibinomial"},
# \code{"scaled_quasibinomial"},
# \code{"ranger"},
# \code{"lm"}.
learn_Q <- function(model_type,
                    history_of_variables,
                    data_learn,
                    formula_strategy = "additive",
                    outcome_name = "weight") {
    max_weight <- max(data_learn[[outcome_name]])
    if (is.null(max_weight) || is.na(max_weight)) {
        stop("The 'weight' column in data_learn must not be NULL or NA.")
    }
    if (max_weight == 0) {
        predict_fun <- function(data) {
            warning("All weights are zero. Returning a constant prediction of zero.")
            rep(0, nrow(data))
        }
        return(predict_fun)
    }
    
    ## String for the formula of the rhs of the formula with the outcome
    if (length(history_of_variables) == 0) {
        history_of_variables_string <- "1"
    } else if (formula_strategy == "additive") {
        history_of_variables_string <- paste(history_of_variables, collapse = "+")
    } else {
        stop("Currently only 'additive' formula strategy is supported.")
    }
    
    if (model_type == "quasibinomial") {
        fit <- stats::glm(
            as.formula(paste0(
                outcome_name,
                " ~ ", history_of_variables_string
            )),
            data = data_learn,
            family = quasibinomial(link = "logit")
        )
        predict_fun <- function(data) {
            predict(fit, data, type = "response")
        }
    } else if (model_type == "scaled_quasibinomial") {
        data_learn$tt_weight <- data_learn[[outcome_name]] / max_weight
        fit <- stats::glm(
            as.formula(paste0(
                "tt_weight ~ ", history_of_variables_string
            )),
            data = data_learn,
            family = quasibinomial
        )
        predict_fun <- function(data) {
            predict(fit, data, type = "response") * max_weight
        }
    } else if (model_type == "ranger") {
        fit <- ranger::ranger(as.formula(paste0(
                           outcome_name,
                           " ~ ", history_of_variables_string
                       )), data = data_learn)
        predict_fun <- function(data) {
            predict(fit, data = data)$predictions
        }  
   } else if (model_type == "lm") {
        fit <- lm(as.formula(paste0(
            outcome_name,
            " ~ ", history_of_variables_string
        )), data = data_learn)
        predict_fun <- function(data) {
            pred <- predict(fit, data, type = "response")
            ## Ensure predictions are non-negative
            pred[pred < 0] <- 0
            return(pred)
        }
   } else if (model_type %in% c("oipcw_expit", "oipcw_probit", "nls_expit", "nls_probit")) {
       Y <- data_learn[[outcome_name]]
       X <- model.matrix(as.formula(paste0(" ~ ", history_of_variables_string)), data = data_learn)

       ## Remove columns of X which are NA in qr.coef(qr(X), Y);
       ## i.e., columns that would normally be removed in a glm fit
       qr_coef <- qr.coef(qr(X),Y)
       X <- X[, !is.na(qr_coef), drop = FALSE]

       ## Start with initial parameters corresponding to intercept only model
       intercept <- ifelse(grepl("expit", model_type), logit(mean(Y)), inv_probit(mean(Y)))
       beta_init <- c(intercept, rep(0, ncol(X) - 1))

       ## Determine if expit or probit from model_type
       link_function <- ifelse(grepl("expit", model_type), expit, probit)

       if (model_type == "oipcw_expit") {
           g <- function(beta, X, Y) {
               eta <- X %*% beta
               as.vector(t(X) %*% (Y - expit(eta)))
           }
       } else if (model_type == "oipcw_probit") {
           g <- function(beta, X, Y) {
               eta <- X %*% beta
               as.vector(t(X) %*% ((Y - probit(eta)) * dnorm(eta) / (probit(eta) * (1 - probit(eta)))))
           }
       } else if (model_type == "nls_expit") {
           g <- function(beta, X, Y) {
               eta <- X %*% beta
               mu <- expit(eta)
               as.vector(t(X) %*% ((Y - mu) * mu * (1 - mu)))
           } 
       } else if (model_type == "nls_probit") {
           g <- function(beta, X, Y) {
               eta <- X %*% beta
               mu <- probit(eta)
               as.vector(t(X) %*% ((Y - mu) * dnorm(eta)))
           }
       }
       
       ## Call to fast C++ implementation of the estimating equation solver (ChatGPT)
       ## only use this with OIPCW;
       ## NLS has serious issues with the C++ solver
       tryCatch({
       if (grepl("oipcw", model_type)) {
           fit <- as.vector(suppressWarnings(estimating_equation_cpp(
               X = X,
               Y = Y,
               model_type = model_type,
               maxit = 1000,
               tol = 1e-8,
               beta = beta_init,
               solve_opts = "force_approx"
           )))
       } else {
           fit <- nleqslv::nleqslv(f = g, x = beta_init, X = X, Y = Y, control = list(maxit = 1000, allowSingular = TRUE))$x
       }},
       error = function(e) {
           warning("The estimating equation solver did not converge: ", e$message)
           fit <<- beta_init
       })
       
       ## Check for NAs or NULL in solution
       if (any(is.na(fit)) || any(is.null(fit))) {
           warning("The estimating equation solver did not converge.")
           fit <- beta_init
       }
      
       ## Check if the solution is very large or if the estimating equation does not seem to be solved
       if (any(abs(fit) > 1e2) ) {
           warning("The solution of the estimating equation solver is very large.")
           message("Estimated parameters / estimating equation value: ")
           print(fit)
           print(g(fit, X, Y))
       }
       if (any(abs(g(fit, X, Y)) > 1e-2)){
           warning("The estimating equation does not seem to be solved which may indicate non-convergence.")
           message("Estimated parameters / estimating equation value: ")
           print(fit)
           print(g(fit, X, Y))
       }

       if (any(is.na(fit))) {
           warning("The estimating equation solver did not converge.")
           fit <- beta_init
       }

       predict_fun <- function(data) {
           X_new <- model.matrix(as.formula(paste0(
               " ~ ", history_of_variables_string
           )), data = data)[, !is.na(qr_coef), drop = FALSE]

           link_function(X_new %*% fit)
       }
   } else if (model_type == "ipcw_glm"){
       stop("ipcw_glm is not implemented yet.")
       ## See Eq. (2) of https://link.springer.com/article/10.1007/s10985-022-09564-6
       ## Can be implemented by fitting a glm with weights;
       ## although we do need censoring survival weights at time $tau$ for that.
   } else {
       ## If flexible, we should pick argmin_(f in cal(F)) sum((Y - g(f(X)))^2, where Y are the outcome weights and g is either expit or probit.
       ## This ensures that predictions will be 0/1-valued.
        formula_w <- paste0(
            outcome_name,
            " ~ ", history_of_variables_string
        )
        predict_fun <- do.call(model_type, list(
                                               character_formula = formula_w,
                                               data = data_learn
                                           ))$predict_fun
    }
    predict_fun
}

## Example function
learn_h2o <- function(character_formula,
                      data,
                      intervened_data = NULL,
                      max_runtime_secs = 30,
                      nfolds = 10,
                      verbose = FALSE,
                      ...) {
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

# coph learner for censoring
learn_coxph <- function(character_formula,
                        data,
                        time_variable = "time"){
  exp_lp <- surv <- hazard <- NULL
  formula_cox <- as.formula(character_formula)
  ## Fit the Cox model
  fit <- coxph(formula_cox, data = data, x = TRUE)
  list(pred = exp(-cumulative_hazard_cox(fit, data, data, time_variable, NULL)$Lambda_minus), fit = fit)
}

learn_glm_logistic <- function(character_formula,
                               data) {
  formula_object <- as.formula(character_formula)
  ## Fit the logistic regression model
  fit <- stats::glm(formula_object, data = data, family = binomial(link = "logit"))
  ## Predict on original data
  list(pred = predict(fit, type = "response"), predict_fun = fit)
}

## Wrapper function to predict the outcome under an intervention
predict_intervention <- function(data, k, predict_fun, static_intervention) {
  event_k <- A_0 <- NULL
  intervened_data <- copy(data)
  if (k > 0) {
    intervened_data[event_k == "A", paste0("A_", k) := static_intervention, env = list(event_k = paste0("event_", k))]
  } else {
    intervened_data[, A_0 := static_intervention]
  }
  f <- predict_fun(intervened_data)
  ## TODO: Check if the predictions are in the range [0,1]
  ## Warn if any predictions are NA or below or above 1
  if (any(is.na(f))) {
    stop("Predictions contain NA values.")
  }
  f
}
