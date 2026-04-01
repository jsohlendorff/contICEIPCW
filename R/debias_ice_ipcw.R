#' @title Computes a one-step estimator of the ICE-IPCW estimator to estimate the mean interventional absolute risk
#' at a given time horizon in continuous time.
#'
#' @param time_horizons A numeric value representing the time horizon at which to estimate the mean interventional absolute risk. 
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
#' @param tmle_update Whether to update the discrete part of the efficient influence function via a TMLE-step instead of one-step.
#' @param penalize_pseudo_outcome Logical; if \code{TRUE}, applies L1-penalization to the regression for the pseudo-outcome. Default is \code{FALSE}.
#' @param penalize_hazard Logical; if \code{TRUE}, applies L1-penalization to the regression for the hazard. Default is \code{FALSE}.
#' @param reduce_colinearity_time Logical; if \code{TRUE}, reduces colinearity in the regression for the pseudo-outcome by using increments of event times. Default is \code{FALSE}.
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
#'  time_horizons = 720,
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
#' time_horizons = 720,
#' model_pseudo_outcome = "lm",
#' model_hazard = "learn_coxph",
#' conservative = TRUE,
#' verbose = TRUE
#' )
#' 
debias_ice_ipcw <- function(
    prepared_data,
    time_horizons = NULL,
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
    tmle_update = FALSE,
    reduce_colinearity_time = FALSE
) {
    if (!inherits(prepared_data, "debiased_prepared")) {
        stop(
            "prepared_data must be of class 'debiased_prepared'. ",
            "Run 'prepare_data' and then 'propensity_scores' first."
        )
    }

    marginal_censoring_fit <- prepared_data$marginal_censoring_fit
    data_info              <- prepared_data$prepared_data_object
    info                   <- data_info$info

    if (is.null(time_horizons)) {
        time_horizons <- info$time_horizon
    }

    data_marginal_censoring <- data_info$data_marginal_censoring
    marginal_censoring      <- data_info$marginal_censoring
    time_covariates         <- data_info$time_covariates
    baseline_covariates     <- data_info$baseline_covariates

    fast_ipcw <- !(
        model_pseudo_outcome %in% c("ipcw_glm_expit", "ipcw_glm_identity")
    )

    if (!fast_ipcw && !marginal_censoring) {
        stop("GLM-based IPCW models require marginal censoring.")
    }

    results <- list()
    ic_list <- list()

    for (v in seq_along(time_horizons)) {
        data <- copy(prepared_data$data)
        th   <- time_horizons[v]

        last_event <- info[time_horizon == th, last_event]
        is_censored <- info[time_horizon == th, is_censored]

        # Remove/rename pooled variables
        pooled_old <- paste0(
            c("time_pooled_", "event_pooled_", "survival_censoring_pooled_"),
            last_event
        )
        pooled_new <- paste0(
            c("time_", "event_", "survival_censoring_"),
            last_event
        )
        set(data, j = paste0(
            c("time_", "event_", "survival_censoring_"),
            last_event
        ), value = NULL)
        setnames(data, old = pooled_old, new = pooled_new)

        set(data,
            j = "ic",
            value = 0
        )
        is_last_event <- TRUE

        # IPCW cumulative weights
        cumulative_inverse_probability_weights(
            data,
            static_intervention,
            th,
            return_ipw,
            last_event
        )

        # Backward iteration over events
        for (k in rev(seq_len(last_event))) {
            set(
                data,
                j = "survival_censoring_k",
                value = data[[paste0("survival_censoring_", k)]]
            )
            set(
                data,
                j = "ipw_cum_weight",
                value = data[[paste0("ipw_cum_weight_", k - 1)]]
            )

            data_at_risk <- get_at_risk_data(data, k, th)
            at_risk_interevent <- data_at_risk$at_risk_interevent
            if (is.null(at_risk_interevent)) next

            data_at_risk <- data_at_risk$at_risk_before_time_horizon

            # Initial q_prediction_prev
            if (is_last_event) {
                data_at_risk[, q_prediction_prev := 0]
            } else {
                data_at_risk <- q_prediction[data_at_risk, on = "id"]
                data_at_risk[is.na(q_prediction_prev), q_prediction_prev := 0]
            }

            # Pseudo-outcomes
            set(
                data_at_risk,
                j = "pseudo_outcome_unweighted",
                value = 1 * (data_at_risk$time_k <= th) * ((data_at_risk$event_k == "Y") + (data_at_risk$event_k %in% c("A", "L")) * data_at_risk$q_prediction_prev)
            )
            set(
                data_at_risk,
                j = "ipcw",
                value = ipcw_k(data_at_risk, k, marginal_censoring_fit, th, is_censored, fast_ipcw, data_at_risk$survival_censoring_k)
            )
            set(
                data_at_risk,
                j = "pseudo_outcome",
                value = data_at_risk$pseudo_outcome_unweighted * data_at_risk$ipcw
            )

            # Regression: q_k
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
                penalize = penalize_pseudo_outcome,
                verbose = verbose,
                reduce_colinearity_time = reduce_colinearity_time,
                time_horizon = th
            )

            # Predict q_k under intervention
            set(
                data_at_risk,
                j = "q_prediction",
                value = predict_intervention(data_at_risk, k-1, q_reg, static_intervention, verbose)
            )

            # Store predictions for next iteration
            q_prediction <- data_at_risk[, .(id, q_prediction_prev = q_prediction)]

            if (anyNA(data_at_risk$q_prediction)) {
                stop("Predictions contain NA values.")
            }

            # Influence-curve updates
            if (!conservative && is_censored) {
                if (tmle_update) {
                    stop("TMLE updates not supported for censored martingales.")
                }
                if (verbose) {
                    message("conservative = FALSE on censored data; proceed with caution.")
                }

                ic_final <- censoring_martingale(
                    data_marginal_censoring,
                    data_at_risk,
                    at_risk_interevent,
                    time_covariates,
                    baseline_covariates,
                    model_hazard,
                    penalize_hazard,
                    model_pseudo_outcome,
                    penalize_pseudo_outcome,
                    th,
                    marginal_censoring,
                    lag,
                    k,
                    grid_size,
                    marginal_censoring_fit,
                    data,
                    static_intervention,
                    reduce_colinearity_time
                )
            } else {
                ic_final <- data[, .(id, ipw_cum_weight)][
                    data_at_risk[, .(id, pseudo_outcome, q_prediction)],
                    on = "id"
                ]

                if (tmle_update) {
                    epsilonhat <- tryCatch(
                        estimating_equation_cpp(
                            X = as.matrix(ic_final$ipw_cum_weight),
                            Y = ic_final$pseudo_outcome,
                            model_type = "oipcw_expit",
                            maxit = 1000,
                            tol = 1e-8,
                            beta = 0,
                            solve_opts = "force_approx",
                            offset = logit(ic_final$q_prediction)
                        )[1, 1],
                        error = function(e) {
                            if (verbose) {
                                message("TMLE update failed: ", e$message)
                            }
                            0
                        }
                    )
                    estimating_equation_tmle_step <- function(epsilon, ipw, pseudo_outcome,q_pred) {
                        as.vector(t(ipw) %*% (pseudo_outcome - expit(logit(q_pred) + epsilonhat * ic_final$ipw_cum_weight)))
                    }
                    g_val <- estimating_equation_tmle_step(epsilonhat, ic_final$ipw_cum_weight, ic_final$pseudo_outcome, ic_final$q_prediction)
                    if (is.na(g_val) || is.null(g_val)) {
                        warning("TMLE update failed to compute estimating equation value; No TMLE update performed for this iteration.")
                        epsilonhat <- 0
                    } else if (abs(g_val) > 1e-1) {
                        warning("TMLE update did not solve the estimating equation with value = ", g_val, "\n. Not updating with TMLE step for this iteration.")
                        epsilonhat <- 0
                    }
                    ## nleqslv::nleqslv(f = g2, x = 0, )
                    ## nleqslv::nleqslv(f = g2, x = 0, ipw = ic_final$ipw_cum_weight, pseudo_outcome = ic_final$pseudo_outcome, q_pred = ic_final$q_prediction, control = list(maxit = 1000, allowSingular = TRUE))$x
                ##  g2(epsilonhat, ic_final$ipw_cum_weight, ic_final$pseudo_outcome, ic_final$q_prediction)
                    ## browser()

                    q_new <- expit(
                        logit(ic_final$q_prediction) +
                        epsilonhat * ic_final$ipw_cum_weight
                    )

                    set(
                        ic_final,
                        j = "q_prediction",
                        value = q_new
                    )
                    set(
                        q_prediction,
                        j = "q_prediction_prev",
                        value = q_new
                    )
                }

                set(
                    ic_final,
                    j = "ipw_cum_weight",
                    value = ic_final$ipw_cum_weight * (ic_final$pseudo_outcome - ic_final$q_prediction)
                )
            }

            ic_final <- ic_final[, .(id, ipw_cum_weight)]

            set(
                data,
                j = "ipw_cum_weight",
                value = NULL
            )
            data <- ic_final[data, on = "id"]
            data[is.na(ipw_cum_weight), ipw_cum_weight := 0]
            set(
                data,
                j = "ic",
                value = data$ic + data$ipw_cum_weight
            )

            is_last_event <- FALSE
        }

        # IPW estimator
        ipw <- if (return_ipw) {
            set(data, j = "ipw", value = 0)
            for (k in seq_len(last_event)) {
                set(
                    data,
                    j = "ipw",
                    value = data[[paste0("ipw_", k)]] + data$ipw
                )
            }
            data[, mean(ipw)]
        } else {
            NA_real_
        }

        # Final estimates
        pred_0 <- predict_intervention(data, 0, q_reg, static_intervention, verbose)
        g_formula_estimate <- data[, mean(pred_0)]

        ic <- data[, ic] + pred_0 - g_formula_estimate

        if (!tmle_update) {
            ice_ipcw_estimate <- g_formula_estimate
            estimate <- ice_ipcw_estimate + mean(ic)
        } else {
            if (verbose) message("Run with tmle_update = FALSE for ICE-IPCW estimate.")
            ice_ipcw_estimate <- NA
            estimate <- g_formula_estimate
        }

        N <- nrow(data)

        results[[v]] <- data.table(
            estimate = estimate,
            se = sd(ic) / sqrt(N),
            lower = estimate - 1.96 * sd(ic) / sqrt(N),
            upper = estimate + 1.96 * sd(ic) / sqrt(N),
            ice_ipcw_estimate = ice_ipcw_estimate,
            ipw = ipw,
            time_horizon = th
        )

        if (return_ic) {
            ic_list[[v]] <- ic
        }
    }

    results <- rbindlist(results)

    if (!return_ic) {
        results
    } else {
        list(results = results, ic = list(ic = ic_list, time_horizons = time_horizons))
    }
}
