### get_propensity_scores.R --- 
#----------------------------------------------------------------------
## Author: Johan Sebastian Ohlendorff
## Created: Feb 26 2026 (17:41) 
## Version: 
## Last-Updated: Mar  6 2026 (13:34) 
##           By: Johan Sebastian Ohlendorff
##     Update #: 221
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
#' @param penalize_treatment Logical; if \code{TRUE}, applies L1 regularization to the treatment propensity score model.
#' @param penalize_hazard Logical; if \code{TRUE}, applies L1 regularization to the hazard model.
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
                              penalize_treatment = FALSE,
                              model_hazard,
                              penalize_hazard = FALSE,
                              lag = NULL,
                              verbose = FALSE) {
    event_number <- id <- ic <- pseudo_outcome <- survival_censoring_k <- event_k <- time_k <- ipw_cum_weight <- ipw_cum_weight_k_prev <- ipw <- ipw_k <- pred_0 <- estimate <- g_formula_estimate <- . <- NULL
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
                               formula_strategy = "additive",
                               penalize = penalize_hazard)
    } else {
        marginal_censoring_fit <- NULL
    }
    ## Unused for now but maybe useful for future extensions
    ## to conservative=FALSE
    ## censoring_models <- list()
    for (k in rev(seq_len(last_event))) {
         ## Create shortcuts for the k'th iteration
         data[, c("event_k", "time_k", "time_k_prev", "event_k_prev", "A_k")
         := list(event_k, time_k, time_k_prev, event_k_prev, A_k),
         env = list(event_k = paste0("event_", k),
                    time_k = paste0("time_", k),
                    event_k_prev = paste0("event_", k - 1),
                    time_k_prev = paste0("time_", k - 1),
                    A_k = paste0("A_", k))]
        
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
                                              time_variable = paste0("time_", k),
                                              penalize = penalize_hazard)

            } else {
                # Ensure Cox.
                # Too difficult without...
                # Have to compute Lambda^c (T_(k,i)- | F_(k-1,i)).
                if (!inherits(marginal_censoring_fit$fit, "coxph")) {
                    stop("Censoring model must be a Cox proportional hazards model when marginal_censoring is TRUE.")
                }

                data_use <- data[event_k_prev %in% c("A", "L")]
                data_use[, c("time", "time_prev") := list(time_k, time_k_prev)]

                data_use <- cumulative_hazard_cox(marginal_censoring_fit$fit, data_use, data_use, time_ref = "time_prev")
                learn_censoring <- list(
                    pred = exp(-data_use$Lambda_minus)
                )
            }
            if (k > 1) {
                data[event_k_prev %in% c("A", "L"),
                     paste0("survival_censoring_", k) := learn_censoring$pred]
            } else {
                data[, paste0("survival_censoring_", k) := learn_censoring$pred]
            }
            ## Currently unused.
            ## censoring_models[[k]] <- learn_censoring$fit
        } else {
            data[, paste0("survival_censoring_", k) := 1]
        }

        ## Fit propensity score (treatment) model
        if (k < last_event) {
            ## check whether all values of A are 1; if so put propensity to 1
            if (all(data[event_k == "A", A_k == 1])) {
                data[event_k == "A", paste0("propensity_",k) := 1]
            } else {
                data[event_k == "A", paste0("propensity_",k) := regression_fit(
                    data = .SD,
                    model_regression = model_treatment,
                    outcome_string = paste0("A_", k),
                    covariates = NULL,
                    formula_strategy = "additive",
                    use_history_of_variables = TRUE,
                    lag = lag,
                    k = k,
                    time_covariates = time_covariates,
                    baseline_covariates = baseline_covariates,
                    type = "propensity",
                    penalize = penalize_treatment
                )]
            }
        }
    }
    ## check whether all values of A_0 are 1; if so put propensity to 1
    if (all(data$A_0 == 1)) {
        data[, propensity_0 := 1]
    } else {
        ## Baseline propensity model
        ## Fit the baseline treatment propensity model
        ## check whethe any baseline covariates should be deleted
        baseline_covariates <- setdiff(
            baseline_covariates,
            names(which(vapply(data[, .SD, .SDcols = baseline_covariates], function(x) length(unique(x)) <= 1, FUN.VALUE = logical(1))))
        )
        data[, propensity_0 := regression_fit(
            data = .SD,
            model_regression = model_treatment,
            outcome_string = "A_0",
            covariates = baseline_covariates,
            formula_strategy = "additive",
            use_history_of_variables = FALSE,
            time_covariates = time_covariates,
            baseline_covariates = baseline_covariates,
            type = "propensity",
            penalize = penalize_treatment
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

