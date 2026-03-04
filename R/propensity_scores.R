### get_propensity_scores.R --- 
#----------------------------------------------------------------------
## Author: Johan Sebastian Ohlendorff
## Created: Feb 26 2026 (17:41) 
## Version: 
## Last-Updated: Mar  4 2026 (22:34) 
##           By: Johan Sebastian Ohlendorff
##     Update #: 185
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
#' @title Function for adding propensity scores (treatment) and censoring models to the prepared data object
#'
#' @param prepared_data An object of class "prepare_data_continuous" containing the prepared data.
#' @param model_treatment A string specifying the type of model to use for the treatment propensity score.
#' Options include \code{"learn_glm_logistic"} (logistic regression).
#' @param model_hazard A string specifying the type of model to use for the cumulative hazard function.
#' Options include \code{"learn_coxph"} (Cox proportional hazards model).
#' @param lag Optional numeric indicating the number of previous events included in the formulas for the models.
#' @param verbose Logical; if \code{TRUE}, prints additional information during model fitting.
#'
#' @export
#' @examples
#' set.seed(15)
#' data_continuous <- simulate_continuous_time_data(
#'   n = 1000,
#'   uncensored = FALSE,
#'   no_competing_events = FALSE,
#'   baseline_rate_list = list(
#'     A = 0.005,
#'     L = 0.001,
#'     C = 0.0008,
#'     Y = 0.0001,
#'     D = 0.00015
#'   )
#' )
#' prep_data <- prepare_data(
#'  data = data_continuous,
#'  max_time_horizon = 720,
#' time_covariates = c("A", "L"),
#' baseline_covariates = c("age", "A_0", "L_0"),
#' marginal_censoring = TRUE
#' )
#' propensity_score_data <- propensity_scores(
#'  prepared_data = prep_data,
#' model_treatment = "learn_glm_logistic",
#' model_hazard = "learn_coxph",
#' verbose = TRUE
#' )
#' 
## Function for getting the propensity scores (treatment) and censoring models
propensity_scores <- function(prepared_data, 
                              model_treatment,
                              model_hazard,
                              lag = NULL,
                              verbose = FALSE) {
    event_number <- id <- ic <- pseudo_outcome <- survival_censoring_k <- event_k <- time_k <- inverse_cumulative_probability_weights <- inverse_cumulative_probability_weights_k_prev <- ipw <- ipw_k <- pred_0 <- estimate <- g_formula_estimate <- . <- NULL
    if (!inherits(prepared_data, "prepare_data_continuous")) {
        stop("prepared_data must be of class 'prepare_data_continuous'.")
    }
    data <- prepared_data$wide_data
    is_censored <- prepared_data$is_censored
    data_marginal_censoring <- prepared_data$data_marginal_censoring
    last_event <- prepared_data$last_event
    marginal_censoring <- prepared_data$marginal_censoring
    time_covariates <- prepared_data$time_covariates
    baseline_covariates <- prepared_data$baseline_covariates
    
    ## Check user input if censored
    if (is_censored && is.null(model_hazard)) {
        stop(
            "Censoring is present, but no censoring model is provided.
             Please provide a censoring model such as `model_hazard = 'learn_coxph'`."
        )
    }

    hazard_minus <- hazard <- event_k <- time_0  <- exp_lp <- surv <- id <- event_k_prev <- survival_censoring_k <- A_k <- propensity_k <- propensity_0 <- NULL
    ## Handle marginal censoring
    if (is_censored && marginal_censoring) {
        ## Remove constant variables
        censoring_covariates <- baseline_covariates[
            data_marginal_censoring[, vapply(.SD, function(x) length(unique(x)) > 1, logical(1)),
                           .SDcols = baseline_covariates]
        ]
        marginal_censoring_fit <- hazard_fit(data = data_marginal_censoring,
                               model_hazard = model_hazard,
                               outcome_string = "Surv(time, event == \"C\")",
                               covariates = censoring_covariates,
                               formula_strategy = "additive")
    } else {
        marginal_censoring_fit <- NULL
    }
    ## Unused for now but maybe useful for future extensions
    ## to conservative=FALSE
    ## censoring_models <- list()
    for (k in rev(seq_len(last_event))) {
        ## Find those at risk of the k'th event; subset data (i.e., people who have not died before the k'th event)
        ## NOTE: For the treatment propensity score, we do not consider
        ## the interarrival times
        data_at_risk <- get_at_risk_data(data, k)
        at_risk_interevent <- data_at_risk$at_risk_interevent
        if (is.null(at_risk_interevent)) {
            next
        }

        ## Fit censoring model if there is censoring
        if (is_censored) {
            if (!marginal_censoring) {
                learn_censoring <- hazard_fit(data = at_risk_interevent,
                                              model_hazard = model_hazard,
                                              outcome_string = paste0("Surv(time_",k,", event_",k," == \"C\")"),
                                              covariates = NULL,
                                              formula_strategy = "additive",
                                              use_history_of_variables = TRUE,
                                              lag = lag,
                                              k = k,
                                              time_covariates = time_covariates,
                                              baseline_covariates = baseline_covariates,
                                              time_variable = paste0("time_", k))

            } else {
                # Ensure Cox.
                # Too difficult without...
                # Have to compute Lambda^c (T_(k,i)- | F_(k-1,i)).
                if (!inherits(marginal_censoring_fit$fit, "coxph")) {
                    stop("Censoring model must be a Cox proportional hazards model when marginal_censoring is TRUE.")
                }

                data_use <- copy(data)[event_prev %in% c("A", "L"), env = list(
                    event_prev = paste0("event_", k - 1)
                )]
                data_use[, c("time", "time_prev") := list(time_k, time_k_prev), env = list(
                    time_k = paste0("time_", k),
                    time_k_prev = paste0("time_", k - 1)
                )]

                data_use <- cumulative_hazard_cox(marginal_censoring_fit$fit, data_use, data_use, time_ref = "time_prev")
                learn_censoring <- list(
                    pred = exp(-data_use$Lambda_minus)
                )
            }
            if (k > 1) {
                data[event_k_prev %in% c("A", "L"),
                     survival_censoring_k := learn_censoring$pred,
                     env = list(
                         survival_censoring_k = paste0("survival_censoring_", k),
                         event_k_prev = paste0("event_", k - 1)
                     )
                     ]
            } else {
                data[, survival_censoring_k := learn_censoring$pred, env = list(
                                                                         survival_censoring_k = paste0("survival_censoring_", k)
                                                                     )]
            }
            ## Currently unused.
            ## censoring_models[[k]] <- learn_censoring$fit
        } else {
            data[, survival_censoring_k := 1, env = list(
                                                  survival_censoring_k = paste0("survival_censoring_", k)
                                              )]
        }

        ## Fit propensity score (treatment) model
        if (k < last_event) {
            history_of_variables_propensity <- get_history_of_variables(
                data[event_k == "A", env = list(event_k = paste0("event_", k))],
                time_covariates,
                baseline_covariates,
                type = "propensity",
                lag = lag,
                k = k
            )
            formula_treatment <- paste0("A_", k, " ~ ", paste0(history_of_variables_propensity, collapse = "+"))
            if (verbose) {
                message("Fitting treatment propensity model for event ", k, " with formula: ", deparse(formula_treatment), "\n")
            }
            ## check whether all values of A are 1; if so put propensity to 1
            if (all(data[event_k == "A", A_k == 1, env = list(
                                                       event_k = paste0("event_", k),
                                                       A_k = paste0("A_", k)
                                                   )])) {
                data[event_k == "A", propensity_k := 1, env = list(
                                                            propensity_k = paste0("propensity_", k),
                                                            event_k = paste0("event_", k)
                                                        )]
            } else {
                data[event_k == "A", propensity_k := withCallingHandlers(
                {
                    do.call(model_treatment, list(
                                                 character_formula = formula_treatment, data = .SD
                                             ))$pred
                },
                error = function(e) {
                    stop("Error in fitting treatment propensity model: ", e, " for event ", k)
                },
                warning = function(w) {
                    message("Warning in fitting treatment propensity model: ", w, " for event ", k)
                }
                ), env = list(
                       propensity_k = paste0("propensity_", k),
                       event_k = paste0("event_", k)
                   )]
            }
        }
    }
    ## check whether all values of A_0 are 1; if so put propensity to 1
    if (all(data$A_0 == 1)) {
        data[, propensity_0 := 1]
    } else {
        ## Baseline propensity model
        formula_treatment <- paste0("A_0 ~ ", paste(
                                                  setdiff(baseline_covariates, "A_0"),
                                                  collapse = "+"
                                              ))
        ## Fit the baseline treatment propensity model
        ## check whethe any baseline covariates should be deleted
        baseline_covariates <- setdiff(
            baseline_covariates,
            names(which(vapply(data[, .SD, .SDcols = baseline_covariates], function(x) length(unique(x)) <= 1, FUN.VALUE = logical(1))))
        )
        if (verbose) {
            message("Fitting baseline treatment propensity model with formula: ", deparse(formula_treatment), "\n")
        }
        data[, propensity_0 := withCallingHandlers(
        {
            do.call(model_treatment, list(
                                         character_formula = formula_treatment, data = .SD
                                     ))$pred
        },
        error = function(e) {
            stop("Error in fitting baseline treatment propensity model: ", e)
        },
        warning = function(w) {
            message("Warning in fitting baseline treatment propensity model: ", w)
        }
        )]
    }
    out<-list(marginal_censoring_fit = marginal_censoring_fit,
              data = data,
              prepared_data_object = prepared_data)
    class(out) <- "debiased_prepared"
    out
    ##list(marginal_censoring_fit = marginal_censoring_fit, censoring_models = censoring_models)
}

######################################################################
### get_propensity_scores.R ends here

