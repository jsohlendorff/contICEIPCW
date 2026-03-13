#' @title Computes a one-step estimator of the ICE-IPCW estimator to estimate the mean interventional absolute risk
#' at a given time horizon in continuous time.
#'
#' @param time_horizon A numeric value representing the time horizon at which to estimate the mean interventional absolute risk. 
#' @param prepared_data An object of class \code{debiased_prepared} containing the prepared data and relevant information for the debiased ICE-IPCW procedure. This object should be obtained by first running the \code{prepare_data} function and then the \code{propensity_scores} function.
#' @param model_pseudo_outcome A string specifying the type of model to use for the iterative conditional expectations estimator.
#'  Options include \code{"tweedie"}, \code{"quasibinomial"}, \code{"scaled_quasibinomial"}, \code{"ranger"}, and \code{"log_normal_mixture"}.
#'  Default is \code{"tweedie"}.
#'  \code{"quasibinomial"} uses a quasi-binomial model. IMPORTANT: Requires outcome between 0 and 1, so cannot be used in the censored case.
#'  \code{"scaled_quasibinomial"} uses a scaled quasi-binomial model, which is similar to \code{"quasibinomial"} but allows for scales outcome to be between 0 and 1
#'  \code{"ranger"} uses a random forest model from the \code{ranger} package.
#'  \code{"log_normal_mixture"} uses a log-normal mixture model, which is useful for continuous outcomes with e.g., allows us to model continuous outcomes with a point mass at 0.
#'
#' @param model_hazard A string specifying the type of model to use for the cumulative hazard function.
#' Options include \code{"learn_coxph"} (Cox proportional hazards model).
#' @param conservative Logical; if \code{TRUE}, do not debias the censoring martingale in the efficient influence function.
#' Results in massive speed up, but slightly less accurate inference.
#' @param return_ipw Logical; if \code{TRUE}, adds inverse probability weight estimator to the output.
#'   Default is \code{TRUE}.
#' @param return_ic Logical; if \code{TRUE}, returns the estimated influence curve (IC) for the ICE-IPCW estimator.
#' @param grid_size Optional numeric indicating the number of grid points to use for `cens_mg_method = "multiple_ice"`.
#' @param lag Optional numeric indicating the number of previous events included in the formulas for the models.
#' @param verbose Logical; if \code{TRUE}, prints additional information during the execution.
#' @param static_intervention Which intervention to consider; either 0 or 1.
#' @param semi_tmle Whether to update the discrete part of the efficient influence function via a TMLE-step instead of one-step.
#' @param penalize_pseudo_outcome Logical; if \code{TRUE}, applies L1-penalization to the regression for the pseudo-outcome. Default is \code{FALSE}.
#' @param penalize_hazard Logical; if \code{TRUE}, applies L1-penalization to the regression for the hazard. Default is \code{FALSE}.
#'
#' @return A named vector containing the following elements:
#' `estimate` - the estimated mean interventional absolute risk at time \code{time_horizon} (debiased)
#' `se` - the standard error of the estimate,
#' `ci_lower` - the lower bound of the confidence interval,
#' `ci_upper` - the upper bound of the confidence interval,
#' `ice_ipcw_estimate` - the ICE-IPCW estimate of the mean interventional absolute risk at time \code{time_horizon} (not debiased),
#' `ipw` - the inverse probability weight estimate (if \code{return_ipw = TRUE}),
#' @details
# ´ Applies inverse probability of censoring weighting (IPCW) to construct an
#' Iterative Conditional Expectation (ICE) estimator for the estimation of a
#' causal effect of a time-varying treatment on a time-to-event outcome
#' that is updated with the efficient influence function, providing inference
#' via a one-step estimator.
#' Current interventions implemented: Static intervention (i.e., intervention at baseline and at each doctor visit to a fixed value).
#' Current target parameters implemented: Absolute risk.
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
#' verbose = TRUE)
#' 
#' result <- debias_ice_ipcw(
#' prepared_data = propensity_score_data,
#' time_horizon = 720,
#' model_pseudo_outcome = "lm",
#' model_hazard = "learn_coxph",
#' conservative = TRUE,
#' verbose = TRUE
#' )
#' 
debias_ice_ipcw <- function(prepared_data,
                            time_horizon,
                            model_pseudo_outcome = "scaled_quasibinomial",
                            penalize_pseudo_outcome = FALSE,
                            model_hazard = "learn_coxph",
                            penalize_hazard = FALSE,
                            conservative = FALSE,
                            static_intervention = 1,
                            return_ipw = TRUE,
                            return_ic = FALSE,
                            grid_size = NULL,
                            lag = NULL,
                            verbose = FALSE,
                            semi_tmle = FALSE) {
    event_number <- id <- ic <- pseudo_outcome <- survival_censoring_k <- event_k <- time_k <- ipw_cum_weight <- ipw_cum_weight_k_prev <- ipw <- ipw_k <- pred_0 <- estimate <- g_formula_estimate <- . <- ipcw <- pseudo_outcome_unweighted <- NULL
    if (!inherits(prepared_data, "debiased_prepared")) {
        stop("prepared_data must be an object of class 'debiased_prepared'.
              Please run the 'prepare_data' function and then the 'propensity_scores'
              function to get an object of class 'debiased_prepared'.")
    }
    marginal_censoring_fit <- prepared_data$marginal_censoring_fit
    data <- prepared_data$data
    data_info <- prepared_data$prepared_data
    is_censored <- data_info$is_censored
    data_marginal_censoring <- data_info$data_marginal_censoring
    last_event <- data_info$last_event
    marginal_censoring <- data_info$marginal_censoring
    time_covariates <- data_info$time_covariates
    baseline_covariates <- data_info$baseline_covariates
    
    data[, ic := 0]    
    is_last_event <- TRUE

    ## IPW weights at each event added to data for the EIF and IPW estimator
    cumulative_inverse_probability_weights(data,
                                           static_intervention,
                                           time_horizon,
                                           return_ipw,
                                           last_event)
    fast_ipcw <- TRUE
    if (model_pseudo_outcome %in% c("ipcw_glm_expit", "ipcw_glm_identity")) {
        fast_ipcw <- FALSE
        ## FIXME
        if (!marginal_censoring) {
            stop("Models with IPCW with glm require marginal censoring to be assumed.")
        }
    }
    
    ## Main procedure for the ICE-IPCW estimator and the one-step update with the efficient influence function
    for (k in rev(seq_len(last_event))) {
        ## Shortcut names
        data[, c("survival_censoring_k", "ipw_cum_weight") :=
        list(survival_censoring_k, ipw_cum_weight_k_prev), env = list(
            survival_censoring_k = paste0("survival_censoring_", k),
            ipw_cum_weight_k_prev = paste0("ipw_cum_weight_", k - 1)
        )]

        data_at_risk <- get_at_risk_data(data, k, time_horizon)
        at_risk_interevent <- data_at_risk$at_risk_interevent
        if (is.null(at_risk_interevent)) {
            next
        }
        data_at_risk <- data_at_risk$at_risk_before_time_horizon
        
        ## Iterated part; use the predictions from the previous iteration
        if (!is_last_event) {
            data_at_risk <- merge(data_at_risk, q_prediction, by = "id", all.x = TRUE)
            data_at_risk[is.na(q_prediction_prev), q_prediction_prev := 0]
        } else {
            data_at_risk[, q_prediction_prev := 0]
        }

        ## Pseudo-outcome tilde(Q)_k
        data_at_risk[, pseudo_outcome_unweighted := 1 * (time_k <= time_horizon) * ((event_k == "Y") + (event_k %in% c("A", "L")) * q_prediction_prev)]
        data_at_risk[, ipcw := ipcw_k(.SD, k, marginal_censoring_fit, time_horizon, is_censored, fast_ipcw, survival_censoring_k)]
        data_at_risk[, pseudo_outcome := pseudo_outcome_unweighted * ipcw]

        ## Fit regression; q_k
        q_reg <- regression_fit(
            data_at_risk,
            model_pseudo_outcome,
            outcome_string = "pseudo_outcome",
            outcome_string_unweighted = "pseudo_outcome_unweighted",
            ipcw_name = "ipcw",
            use_history_of_variables = TRUE,
            lag = lag,
            k = k,
            time_covariates = time_covariates,
            baseline_covariates = baseline_covariates,
            type = "pseudo_outcome",
            penalize = penalize_pseudo_outcome
        )
        
        ## Predict q_k under the previous intervention
        data_at_risk[, q_prediction := predict_intervention(.SD, k-1, q_reg, static_intervention)]

        ## Save values for next iteration
        q_prediction <- data_at_risk[, c("q_prediction", "id"), with = FALSE]
        setnames(q_prediction, "q_prediction", "q_prediction_prev")

        ## Throw error if any predictions are NA
        if (any(is.na(data_at_risk$q_prediction))) {
            stop("Predictions contain NA values.")
        }

        if (!conservative & is_censored) {
            if (semi_tmle) stop("semi-tmle not implemented yet for censored martingale")
            message("conservative = FALSE on censored data. You're on shaky ground...")
            ic_final <- censoring_martingale(data_marginal_censoring,
                                             data_at_risk,
                                             at_risk_interevent,
                                             time_covariates,
                                             baseline_covariates,
                                             model_hazard,
                                             penalize_hazard,
                                             model_pseudo_outcome,
                                             penalize_pseudo_outcome,
                                             time_horizon,
                                             marginal_censoring,
                                             lag,
                                             k,
                                             grid_size,
                                             marginal_censoring_fit,
                                             data,
                                             static_intervention)
        } else {
            ## If conservative, we do not compute the martingale terms
            ic_final <- merge(data_at_risk[, c("pseudo_outcome", "q_prediction", "id")], data[, c("ipw_cum_weight", "id")], by = "id")
            if (semi_tmle) {
                ## Note: Solving the equation for scaled q_predictionictions and scaled pseudo_outcomes, correspond to getting epsilon from original problem
                tryCatch({
                    epsilonhat <- estimating_equation_cpp(
                        X = as.matrix(ic_final$ipw_cum_weight),
                        Y = ic_final$pseudo_outcome,
                        model_type = "oipcw_expit",
                        maxit = 1000,
                        tol = 1e-8,
                        beta = 0,
                        solve_opts = "force_approx",
                        offset = logit(ic_final$q_prediction)
                    )[1,1]
                }, error = function(e) {
                    warning("Error in glm.fit for TMLE update. Setting epsilonhat to 0.")
                    epsilonhat <<- 0
                })
                ##  # Debug
                ##  g2 <- function(epsilon, ipw, pseudo_outcome,q_pred) {
                ##  as.vector(t(ipw) %*% (pseudo_outcome - expit(logit(q_pred) + epsilonhat * ic_final$ipw_cum_weight)))
                ##  }
                ##  g2(epsilonhat, ic_final$ipw_cum_weight, ic_final$pseudo_outcome, ic_final$q_prediction)
                q_prediction_prev <- expit(logit(ic_final$q_prediction) + epsilonhat * (ic_final$ipw_cum_weight))
                ic_final$q_prediction <- q_prediction_prev
                q_prediction$q_prediction_prev <- q_prediction_prev
            }
            ic_final <- ic_final[, ipw_cum_weight := ipw_cum_weight * (pseudo_outcome - q_prediction)]
        }
        ic_final <- ic_final[, c("ipw_cum_weight", "id")]
        ## Now add the influence curve to the data 
        data[, ipw_cum_weight := NULL]
        data <- merge(ic_final, data, by = "id", all = TRUE)
        data[is.na(ipw_cum_weight), ipw_cum_weight := 0]

        data[, ic := ic + ipw_cum_weight]
        is_last_event <- FALSE
    }
    if (return_ipw) {
        data[, ipw := 0]
        for (k in seq_len(last_event)) {
            data[, ipw := ipw + ipw_k, env = list(ipw_k = paste0("ipw_", k))]
        }
        data[, ipw := mean(ipw)]
    }
    data[, pred_0 := predict_intervention(.SD, 0, q_reg, static_intervention)]
    data[, g_formula_estimate := mean(pred_0)]
    data[, ic := ic + pred_0 - g_formula_estimate]
    if (!semi_tmle) {
        data[, estimate := g_formula_estimate + mean(ic)]
        result <- data[, .(
            estimate = estimate[.N],
            se = sd(ic) / sqrt(.N),
            lower = estimate[.N] - 1.96 * sd(ic) / sqrt(.N),
            upper = estimate[.N] + 1.96 * sd(ic) / sqrt(.N),
            ice_ipcw_estimate = g_formula_estimate[.N],
            ipw = ipw[.N]
        )]
    } else {
        message("Rerun with semi_tmle = FALSE to get 'ice_ipcw_estimate'")
        result <- data[, .(
            estimate = g_formula_estimate[.N],
            se = sd(ic) / sqrt(.N),
            lower = g_formula_estimate[.N] - 1.96 * sd(ic) / sqrt(.N),
            upper = g_formula_estimate[.N] + 1.96 * sd(ic) / sqrt(.N),
            ice_ipcw_estimate = NA,
            ipw = ipw[.N])]
    }
    if (return_ic) {
        list(result = result, ic = data[, ic])
    } else {
        result
    }
}
## TODO: Add possibility to use IPW as the last regression when few event points are available
## TODO: Add possibility to simulate (impute) when few event points are available
## TODO: Add cross-fitting as a possibility
## TODO: Add pooling later
