### propensity_scores.R --- 
#----------------------------------------------------------------------
## Author: Johan Sebastian Ohlendorff
## Created: Feb 26 2026 (17:41) 
## Version: 
## Last-Updated: May 13 2026 (12:32) 
##           By: Johan Sebastian Ohlendorff
##     Update #: 487
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
#' @param static_intervention Numeric value indicating the treatment level for the static intervention (default is 1).
#' @param exclude_latest_covariate Optional character vector of covariate names to exclude the latest value of in the propensity score models.
#' @param reduce_colinearity_time Logical; if \code{TRUE}, reduces colinearity in the regression for treatment by using increments of event times. Default is \code{FALSE}.
#' @param gbound Optional numeric value for lower bounding the propensity scores to avoid extreme weights. If \code{NULL}, no bounding is applied.
#' @param save_coefficients Logical; if \code{TRUE}, saves the coefficients of the fitted models in the output object. Default is \code{FALSE}.
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
                              static_intervention = 1,
                              exclude_latest_covariate = NULL,
                              reduce_colinearity_time = FALSE,
                              gbound = NULL,
                              verbose = FALSE,
                              save_coefficients = FALSE) {
    if (!inherits(prepared_data, "prepare_data_continuous")) {
        stop("prepared_data must be of class 'prepare_data_continuous'.")
    }

    data <- prepared_data$wide_data    
    data_marginal_censoring <- prepared_data$data_marginal_censoring
    info <- prepared_data$info
    time_covariates <- prepared_data$time_covariates
    baseline_covariates <- prepared_data$baseline_covariates
    marginal_censoring <- prepared_data$marginal_censoring

    ## Check user input if censored
    if (any(info$is_censored) && is.null(model_hazard)) {
        stop("Censoring is present, but no censoring model is provided.")
    }

    ## Handle marginal censoring
    if (any(info$is_censored) && marginal_censoring) {

        ## Remove constant variables
        constant_vars <- vapply(
            data_marginal_censoring[, ..baseline_covariates],
            function(x) length(unique(x)) > 1,
            logical(1)
        )
        censoring_covariates <- baseline_covariates[constant_vars]

        marginal_censoring_fit <- hazard_fit(
            data = data_marginal_censoring,
            model_hazard = model_hazard,
            outcome_string = "Surv(time, event == \"C\")",
            covariates = censoring_covariates,
            formula_strategy = "additive",
            penalize = penalize_hazard,
            verbose = verbose
        )
    } else {
        marginal_censoring_fit <- NULL
    }

    is_censored_max <- any(info$is_censored)
    last_event <- max(info$last_event)
    unique_last_events <- unique(info$last_event)

    if (save_coefficients) {
        coefficients_list <- list()
        if (marginal_censoring) {
            coefficients_list$marginal_censoring <- coef(marginal_censoring_fit$fit)
        }
        if (verbose) {
            message("Coefficients for non-marginal censoring models will not be saved")
        }
        coefficients_propensity <- list()
    }

    ############################################################################
    ## MAIN LOOP OVER EVENTS k
    ############################################################################
    for (k in rev(seq_len(last_event))) {        
        cols <- c("event_k", "time_k", "time_k_prev", "event_k_prev", "A_k")
        vals <- list(
            data[[paste0("event_", k)]],
            data[[paste0("time_", k)]],
            data[[paste0("time_", k - 1)]],
            data[[paste0("event_", k - 1)]],
            data[[paste0("A_", k)]]
        )

        for (j in seq_along(cols)) {
            set(data, j = cols[j], value = vals[[j]])
        }

        ## Determine risk set
        data_at_risk <- get_at_risk_data(data, k)
        at_risk_interevent <- data_at_risk$at_risk_interevent
        if (is.null(at_risk_interevent)) next

        ########################################################################
        ## CENSORING MODEL
        ########################################################################
        if (is_censored_max) {

            censoring_wrapper <- function(data2,
                                          at_risk_interevent,
                                          model_hazard,
                                          outcome_string,
                                          covariates,
                                          formula_strategy,
                                          use_history_of_variables,
                                          lag,
                                          k,
                                          time_covariates,
                                          baseline_covariates,
                                          penalize,
                                          verbose,
                                          time_variable,
                                          event_variable) {

                if (!marginal_censoring) {

                    learn_censoring <- hazard_fit(
                        data = at_risk_interevent,
                        model_hazard = model_hazard,
                        outcome_string = paste0("Surv(", time_variable, "_", k,
                                                ", ", event_variable, "_", k,
                                                " == \"C\")"),
                        covariates = NULL,
                        formula_strategy = "additive",
                        use_history_of_variables = TRUE,
                        lag = lag,
                        k = k,
                        time_covariates = time_covariates,
                        baseline_covariates = baseline_covariates,
                        time_variable = paste0(time_variable, "_", k),
                        penalize = penalize_hazard,
                        verbose = verbose
                    )
                    learn_censoring$pred

                } else {

                    if (!inherits(marginal_censoring_fit$fit, "coxph")) {
                        stop("Censoring model must be Cox when marginal_censoring=TRUE.")
                    }

                    data_use <- data2[event_k_prev %in% c("A","L")]

                    set(data_use, j = "time",
                        value = data_use[[paste0(time_variable, "_", k)]])
                    set(data_use, j = "time_prev",
                        value = data_use[[paste0("time_", k - 1)]])

                    data_use2 <- cumulative_hazard_cox(
                        fit = marginal_censoring_fit$fit,
                        data = data_use,
                        covariate_data = data_use,
                        time_ref = "time_prev",
                        baseline_hazard = marginal_censoring_fit$baseline_hazard
                    )

                    pred <- exp(-data_use2$Lambda_minus)

                    if (any(is.na(pred))) {
                        stop("NA values for IPCW.")
                    }
                    if (any(pred == 0)) {
                        stop("Zero values for IPCW.")
                    }

                    pred
                }
            }

            pred <- censoring_wrapper(
                at_risk_interevent = at_risk_interevent,
                data2 = data,
                model_hazard = model_hazard,
                outcome_string = NULL,
                covariates = NULL,
                formula_strategy = "additive",
                use_history_of_variables = TRUE,
                lag = lag,
                k = k,
                time_covariates = time_covariates,
                baseline_covariates = baseline_covariates,
                penalize = penalize_hazard,
                verbose = verbose,
                time_variable = "time",
                event_variable = "event"
            )

            colname <- paste0("survival_censoring_", k)

            if (k > 1) {
                rows <- which(data$event_k_prev %in% c("A","L"))
                set(data, i = rows, j = colname, value = pred)
            } else {
                set(data, j = colname, value = pred)
            }

            ## Pooled censoring for unique last events
            if (k %in% unique_last_events) {
                pooled_col <- paste0("survival_censoring_pooled_", k)
                if (k < last_event) {
                    pred2 <- censoring_wrapper(
                        at_risk_interevent = at_risk_interevent,
                        data2 = data,
                        model_hazard = model_hazard,
                        outcome_string = NULL,
                        covariates = NULL,
                        formula_strategy = "additive",
                        use_history_of_variables = TRUE,
                        lag = lag,
                        k = k,
                        time_covariates = time_covariates,
                        baseline_covariates = baseline_covariates,
                        penalize = penalize_hazard,
                        verbose = verbose,
                        time_variable = "time_pooled",
                        event_variable = "event_pooled"
                    )

                    if (k > 1) {
                        rows <- which(data$event_k_prev %in% c("A","L"))
                        set(data, i = rows, j = pooled_col, value = pred2)
                    } else {
                        set(data, j = pooled_col, value = pred2)
                    }
                } else {
                    basecol <- paste0("survival_censoring_", k)
                    set(data, j = pooled_col, value = data[[basecol]])
                }
            }

        } else {
            if (k %in% unique_last_events) {
                set(data, j = paste0("survival_censoring_pooled_", k), value = 1)
            }
            set(data, j = paste0("survival_censoring_", k), value = 1)
        }

        ########################################################################
        ## PROPENSITY MODEL FOR A_k (treatment)
        ########################################################################
        if (k < last_event) {

            A_k_col <- paste0("A_", k)
            pcol <- paste0("propensity_", k)
            rowsA <- which(data$event_k == "A")

            if (all(data[rowsA][[A_k_col]] == 1)) {

                set(data, i = rowsA, j = pcol, value = 1)

            } else {

                varcol <- paste0("A_", k, "_var")
                set(data, j = varcol,
                    value = data[[A_k_col]] == static_intervention)

                preds <- regression_fit(
                    data = data[rowsA],
                    model_regression = model_treatment,
                    outcome_string = varcol,
                    covariates = NULL,
                    formula_strategy = "additive",
                    use_history_of_variables = TRUE,
                    lag = lag,
                    k = k,
                    time_covariates = time_covariates,
                    baseline_covariates = baseline_covariates,
                    type = "propensity",
                    penalize = penalize_treatment,
                    exclude_latest_covariate = exclude_latest_covariate,
                    verbose = verbose,
                    reduce_colinearity_time = reduce_colinearity_time,
                    save_coefficients = save_coefficients
                )
                if (save_coefficients) {
                    coefficients_propensity[[paste0("A_", k)]] <- preds$coefficients
                }
                preds <- preds$predictions
                set(data, j = varcol, value = NULL)
                if (any(is.na(preds))) {
                    stop("NA in propensity scores for event ", k)
                }
                if (!is.null(gbound)) {
                    set(data, i = rowsA, j = pcol, value = pmax(data[rowsA][[pcol]], gbound))
                } else {
                    set(data, i = rowsA, j = pcol, value = preds)
                }
            }
        }
    }  ## end loop over k

    ############################################################################
    ## BASELINE TREATMENT MODEL A_0
    ############################################################################
    if (all(data$A_0 == 1)) {

        set(data, j = "propensity_0", value = 1)

    } else {
        baseline_covariates <- setdiff(baseline_covariates, "A_0")
        keep <- !vapply(
            data[, ..baseline_covariates],
            function(x) length(unique(x)) <= 1,
            logical(1)
        )
        baseline_covariates2 <- baseline_covariates[keep]

        set(data, j = "A_0_var", value = data$A_0 == static_intervention)

        preds0 <- regression_fit(
            data = data,
            model_regression = model_treatment,
            outcome_string = "A_0_var",
            covariates = baseline_covariates2,
            formula_strategy = "additive",
            use_history_of_variables = FALSE,
            time_covariates = time_covariates,
            baseline_covariates = baseline_covariates2,
            type = "propensity",
            penalize = penalize_treatment,
            exclude_latest_covariate = exclude_latest_covariate,
            verbose = verbose,
            save_coefficients = save_coefficients
        )
        if (save_coefficients) {
            coefficients_propensity[["A_0"]] <- preds0$coefficients
        }
        preds0 <- preds0$predictions

        set(data, j = "A_0_var", value = NULL)

        if (any(is.na(preds0))) {
            stop("NA values in baseline propensity scores.")
        }
        
        if (!is.null(gbound)) {
            set(data, j = "propensity_0", value = pmax(preds0, gbound))
        } else {
            set(data, j = "propensity_0", value = preds0)
        }
    }

    if (save_coefficients) {
        coefficients_list[["propensity"]] <- coefficients_propensity
    } else {
        coefficients_list <- NULL
    }

    ############################################################################
    ## Remove temporary columns
    ############################################################################
    temp_cols <- c("event_k", "time_k", "time_k_prev", "event_k_prev", "A_k")
    for (cc in temp_cols) set(data, j = cc, value = NULL)

    out <- list(
        marginal_censoring_fit = marginal_censoring_fit,
        data = data,
        prepared_data_object = prepared_data,
        coefficients_propensity = coefficients_list
    )
    class(out) <- "debiased_prepared"
    out
}
######################################################################
### get_propensity_scores.R ends here

