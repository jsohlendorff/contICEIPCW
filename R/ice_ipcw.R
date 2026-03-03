#' @title Computes a one-step estimator of the ICE-IPCW estimator to estimate the mean interventional absolute risk
#' at a given time horizon in continuous time.
#'
#' @param data An object containing two data frames: \code{baseline_data} and
#'  \code{timevarying_data}. The \code{baseline_data} should contain baseline
#'  covariates, that is `id`, `A_0`, `L_0`, and the treatment variable
#'  `A_0` must be a binary variable (0/1) and `L_0`
#'  (the initial value of the time-varying covariates) can only be one-dimesional.
#'  Additional baseline covariates that are not time-varying can be added here.
#'  The \code{timevarying_data} should contain
#'  `id`, `time`, `event`, `A`, and `L` columns, where `event` is a
#'  factor with levels (`A`, `L`, `C`, `Y`, `D`), i.e., which event happened.
#'  `time` is the time of the corresponding event, `A` is the time-varying treatment variable
#'  and `L` is the time-varying covariate which must be one-dimesional
#'  Note that `A`, `L`, `C`, `Y`, and `D` are the event types, corresponding to
#'  `Y` (event of interest), `D` (competing event), `A` (visiting event), `L` (covariate event),
#'  `C` (censoring event).
#' @param time_horizon A numeric value representing the time horizon for the analysis.
#' @param model_pseudo_outcome A string specifying the type of model to use for the iterative conditional expectations estimator.
#'  Options include \code{"tweedie"}, \code{"quasibinomial"}, \code{"scaled_quasibinomial"}, \code{"ranger"}, and \code{"log_normal_mixture"}.
#'  Default is \code{"tweedie"}.
#'  \code{"quasibinomial"} uses a quasi-binomial model. IMPORTANT: Requires outcome between 0 and 1, so cannot be used in the censored case.
#'  \code{"scaled_quasibinomial"} uses a scaled quasi-binomial model, which is similar to \code{"quasibinomial"} but allows for scales outcome to be between 0 and 1
#'  \code{"ranger"} uses a random forest model from the \code{ranger} package.
#'  \code{"log_normal_mixture"} uses a log-normal mixture model, which is useful for continuous outcomes with e.g., allows us to model continuous outcomes with a point mass at 0.
#'
#' @param model_treatment A string specifying the type of model to use for the treatment propensity score.
#' Options include \code{"learn_glm_logistic"} (logistic regression).
#' @param model_hazard A string specifying the type of model to use for the cumulative hazard function.
#' Options include \code{"learn_coxph"} (Cox proportional hazards model).
#' @param conservative Logical; if \code{TRUE}, do not debias the censoring martingale in the efficient influence function.
#' Results in massive speed up, but slightly less accurate inference.
#' @param time_covariates A character vector of column names in \code{data} that are
#'   treated as time-varying covariates. Must include values of time-varying covariates at baseline.
#' @param baseline_covariates A character vector of column names in \code{data} that are
#'   considered baseline (time-invariant) covariates. Must include treatment and time-varying covariates.
#' @param last_event Optional numeric indicating the last nonterminal event number to consider
#'   in the outcome.
#' @param return_ipw Logical; if \code{TRUE}, adds inverse probability weight estimator to the output.
#'   Default is \code{TRUE}.
#' @param return_ic Logical; if \code{TRUE}, returns the estimated influence curve (IC) for the ICE-IPCW estimator.
#' @param grid_size Optional numeric indicating the number of grid points to use for `cens_mg_method = "multiple_ice"`.
#' @param lag Optional numeric indicating the number of previous events included in the formulas for the models.
#' @param verbose Logical; if \code{TRUE}, prints additional information during the execution.
#' @param marginal_censoring Logical; if \code{TRUE}, assumes censoring depends only on baseline covariates.
#' @param static_intervention Which intervention to consider; either 0 or 1.
#' @param semi_tmle Whether to update the discrete part of the efficient influence function via a TMLE-step instead of one-step.
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
#'
#' debias_ice_ipcw(
#'   data = data_continuous,
#'   time_horizon = 720,
#'   model_pseudo_outcome = "oipcw_expit",
#'   model_treatment = "learn_glm_logistic",
#'   model_hazard = "learn_coxph",
#'   time_covariates = c("A", "L"),
#'   baseline_covariates = c("age", "A_0", "L_0"),
#'   conservative = TRUE
#' )
debias_ice_ipcw <- function(data,
                            time_horizon,
                            model_pseudo_outcome = "scaled_quasibinomial",
                            model_treatment = "learn_glm_logistic",
                            model_hazard = "learn_coxph",
                            marginal_censoring = FALSE,
                            conservative = FALSE,
                            time_covariates,
                            baseline_covariates,
                            last_event = NULL,
                            static_intervention = 1,
                            return_ipw = TRUE,
                            return_ic = FALSE,
                            grid_size = NULL,
                            lag = NULL,
                            verbose = FALSE,
                            semi_tmle = FALSE) {
    event_number <- id <- ic <- pseudo_outcome <- survival_censoring_k <- event_k <- time_k <- inverse_cumulative_probability_weights <- inverse_cumulative_probability_weights_k_prev <- ipw <- ipw_k <- pred_0 <- estimate <- g_formula_estimate <- . <- NULL
    ## Check user input
    check_input(baseline_covariates, time_covariates, data, time_horizon)

    ## Get timevarying data and baseline data and add event number by id
    timevarying_data <- data$timevarying_data[, event_number := seq_len(.N), by = id]
    baseline_data <- data$baseline_data

    ## If last event number not provided,
    ## select last event number adaptively because the iterative
    ## regressions may not have sufficient data to fit the models for later events.
    ## NOTE: Modifies data.
    last_event <- select_last_event(timevarying_data, time_horizon, last_event)

    censoring_info_result <- censoring_info(timevarying_data,
                                            baseline_data,
                                            time_horizon,
                                            marginal_censoring,
                                            model_hazard,
                                            model_pseudo_outcome)
    is_censored <- censoring_info_result$is_censored
    data_baseline <- censoring_info_result$data_baseline
    
    ## Convert the data from long format to wide format
    data <- widen_continuous_data(timevarying_data, baseline_data, time_covariates)
    data[, ic := 0]    
    is_last_event <- TRUE

    ## Get propensity scores and models for the censoring.
    ## NOTE: Modifies data in place, so that the propensity scores are added to the data.
    marginal_censoring_fit <- tryCatch(
    {
        propensity_scores(
            last_event,
            data,
            time_horizon,
            model_treatment,
            model_hazard,
            is_censored,
            time_covariates,
            baseline_covariates,
            marginal_censoring,
            lag,
            verbose,
            data_baseline
        )
    },
    error = function(e) {
        stop("Error in getting censoring/propensity models: ", e)
    })

    ## IPW weights at each event added to data for the EIF and IPW estimator
    cumulative_inverse_probability_weights(data,
                                           static_intervention,
                                           time_horizon,
                                           return_ipw,
                                           last_event)
    
    ## Main procedure for the ICE-IPCW estimator and the one-step update with the efficient influence function
    for (k in rev(seq_len(last_event))) {
        data_at_risk <- get_at_risk_data(data, k, time_horizon)
        at_risk_interevent <- data_at_risk$at_risk_interevent
        if (is.null(at_risk_interevent)) {
            next
        }
        data_at_risk <- data_at_risk$at_risk_before_time_horizon
        
        ## Iterated part; use the predictions from the previous iteration
        if (!is_last_event) {
            data_at_risk <- merge(data_at_risk, q_prediction, by = "id", all.x = TRUE)
            data_at_risk[is.na(future_prediction), future_prediction := 0]            
        } else {
            data_at_risk[, future_prediction := 0]
        }

        ## Pseudo-outcome tilde(Q)_k
        data_at_risk[, pseudo_outcome := 1 / (survival_censoring_k) * ((event_k == "Y" & time_k <= time_horizon) + (event_k %in% c("A", "L")) * future_prediction), env = list(
                                                                                                                                                             survival_censoring_k = paste0("survival_censoring_", k),
                                                                                                                                                             event_k = paste0("event_", k),
                                                                                                                                                             time_k = paste0("time_", k)
                                                                                                                                                             )]

        ## Fit regression; q_k
        q_reg <- regression_fit(
            data_at_risk,
            model_pseudo_outcome,
            outcome_string = "pseudo_outcome",
            use_history_of_variables = TRUE,
            lag = lag,
            k = k,
            time_covariates = time_covariates,
            baseline_covariates = baseline_covariates,
            type = "pseudo_outcome"
        )
        
        ## Predict q_k under the previous intervention
        data_at_risk[, q_prediction := predict_intervention(.SD, k-1, q_reg, static_intervention)]

        ## Save values for next iteration
        q_prediction <- data_at_risk[, c("q_prediction", "id"), with = FALSE]
        setnames(q_prediction, "q_prediction", "future_prediction")

        ## Throw error if any predictions are NA
        if (any(is.na(data_at_risk$q_prediction))) {
            stop("Predictions contain NA values.")
        }
        data_at_risk[, inverse_cumulative_probability_weights := inverse_cumulative_probability_weights_k_prev, env = list(inverse_cumulative_probability_weights_k_prev = paste0("inverse_cumulative_probability_weights_", k - 1))]
        data[, inverse_cumulative_probability_weights := inverse_cumulative_probability_weights_k_prev, env = list(inverse_cumulative_probability_weights_k_prev = paste0("inverse_cumulative_probability_weights_", k - 1))]

        if (!conservative & is_censored) {
            if (semi_tmle) stop("semi-tmle not implemented yet for censored martingale")
            ic_final <- censoring_martingale(data_baseline,
                                             data_at_risk,
                                             at_risk_interevent,
                                             time_covariates,
                                             baseline_covariates,
                                             model_hazard,
                                             model_pseudo_outcome,
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
            ic_final <- merge(data_at_risk[, c("pseudo_outcome", "q_prediction", "id")], data[, c("inverse_cumulative_probability_weights", "id")], by = "id")
            if (semi_tmle) {
                max_pseudo_outcome <- max(ic_final$pseudo_outcome)
                ## Note: Solving the equation for scaled q_predictionictions and scaled pseudo_outcomes, correspond to getting epsilon from original problem
                ic_final$f_pseudo_outcome <- ic_final$pseudo_outcome / max_pseudo_outcome
                ic_final$f_q_prediction <- ic_final$q_prediction / max_pseudo_outcome
                epsilonhat <- stats::glm(f_pseudo_outcome~inverse_cumulative_probability_weights-1+offset(qlogis(f_q_prediction)),family="quasibinomial",data = ic_final)$coefficients[1]
                ic_final[, c("f_pseudo_outcome", "f_q_prediction") := NULL]
                future_prediction <- stats::plogis(stats::qlogis(ic_final$q_prediction) + epsilonhat * (ic_final$inverse_cumulative_probability_weights))
                ic_final$q_prediction <- future_prediction
                q_prediction$future_prediction <- future_prediction
            }
            ic_final <- ic_final[, inverse_cumulative_probability_weights := inverse_cumulative_probability_weights * (pseudo_outcome - q_prediction)] # pseudo_outcome: Z. pred: Q
        }
        ic_final <- ic_final[, c("inverse_cumulative_probability_weights", "id")]
        ## Now add the influence curve to the data data
        data[, inverse_cumulative_probability_weights := NULL]
        data <- merge(ic_final, data, by = "id", all = TRUE)
        data[is.na(inverse_cumulative_probability_weights), inverse_cumulative_probability_weights := 0]

        data[, ic := ic + inverse_cumulative_probability_weights]
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
        result <- data[, .(
            estimate = g_formula_estimate[.N],
            se = sd(ic) / sqrt(.N),
            lower = g_formula_estimate[.N] - 1.96 * sd(ic) / sqrt(.N),
            upper = g_formula_estimate[.N] + 1.96 * sd(ic) / sqrt(.N),
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
