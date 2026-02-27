### get_propensity_scores.R --- 
#----------------------------------------------------------------------
## Author: Johan Sebastian Ohlendorff
## Created: Feb 26 2026 (17:41) 
## Version: 
## Last-Updated: Feb 27 2026 (19:25) 
##           By: Johan Sebastian Ohlendorff
##     Update #: 42
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
## Function for getting the propensity scores (treatment) and censoring models
propensity_scores <- function(last_event,
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
                              data_baseline) {
    ## Handle marginal censoring
    if (is_censored && marginal_censoring) {
        ## Remove constant variables
        censoring_covariates <- baseline_covariates[
            data_baseline[, vapply(.SD, function(x) length(unique(x)) > 1, logical(1)),
                           .SDcols = baseline_covariates]
        ]
        marginal_censoring_fit <- hazard_fit(data = data_baseline,
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
                                              baseline_covariates = baseline_covariates)

            } else {
                ## FIXME: make marginal censoring hazard work without cox
                if (!inherits(marginal_censoring_fit$fit, "coxph")) {
                    stop("Censoring model must be a Cox proportional hazards model when marginal_censoring is TRUE.")
                }
                base_hazard <- as.data.table(basehaz(marginal_censoring_fit$fit, centered = FALSE))
                ## First find survival probability of not being censored at time T_k- (when the covariates are zero)
                baseline_hazard_minus <- copy(base_hazard)
                baseline_hazard_minus$hazard <- c(0, head(baseline_hazard_minus$hazard, -1))
                baseline_hazard <- copy(base_hazard)
                setnames(baseline_hazard_minus, c("time", "hazard"), c(paste0("time_", k), "hazard_minus"))
                setnames(baseline_hazard, "time", paste0("time_", k))
                baseline_hazard_minus <- baseline_hazard_minus[data[, c("id", paste0("event_", k), paste0("time_", k), baseline_covariates), with = FALSE], roll = TRUE, on = paste0("time_", k)]
                baseline_hazard_minus[is.na(hazard_minus), hazard_minus := 0]
                baseline_hazard <- baseline_hazard[data[, c("id", paste0("event_", k), paste0("time_", k), baseline_covariates), with = FALSE], roll = TRUE, on = paste0("time_", k)]
                baseline_hazard[is.na(hazard), hazard := 0]
                baseline_hazard_minus <- merge(baseline_hazard, baseline_hazard_minus, by = c("id", baseline_covariates, paste0("event_", k), paste0("time_", k)))
                ## Split into whether or not, the event time is registered on the data set consisting only of terminal events
                ## For the Cox model, if t is not an event time, then the survival fucntion G(t- | x) = G(t | x)
                ## A and L events are NOT event times for the marginal Cox model
                baseline_hazard_minus[, hazard_minus := hazard_minus * 1 * (event_k %in% c("C", "Y", "D", "tauend")) + hazard * 1 * (event_k %in% c("A", "L")), env = list(event_k = paste0("event_", k))]
                baseline_hazard_minus[, c(paste0("event_", k), "hazard") := NULL]

                ## Then find survival probability of not being censored at time T_(k-1) (when the covariates are zero)
                if (k > 1) {
                    baseline_hazard <- copy(base_hazard)
                    setnames(baseline_hazard, "time", paste0("time_", k - 1))
                    baseline_hazard <- baseline_hazard[data[, c("id", paste0("time_", k - 1), baseline_covariates), with = FALSE], roll = TRUE, on = paste0("time_", k - 1)]
                    baseline_hazard[is.na(hazard), hazard := 0]
                } else {
                    baseline_hazard <- baseline_hazard_minus[, -c("hazard_minus", paste0("time_", k)), with = FALSE]
                    baseline_hazard[, time_0 := 0]
                    baseline_hazard[, hazard := 0]
                }
                baseline_hazard <- merge(baseline_hazard, baseline_hazard_minus, by = c("id", baseline_covariates))
                baseline_hazard[, exp_lp := predict(marginal_censoring_fit$fit, newdata = .SD, type = "risk", reference = "zero")]

                ## Survival probability of not being censored at time T_k- (when the covariates are zero) given not censored at time T_(k-1)
                baseline_hazard[, surv := exp(-exp_lp * (hazard_minus - hazard))]
                learn_censoring <- list()
                learn_censoring$pred <- baseline_hazard[id %in% at_risk_interevent$id, surv]
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
            names(which(sapply(data[, ..baseline_covariates], function(x) length(unique(x)) <= 1)))
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
    marginal_censoring_fit
    ##list(marginal_censoring_fit = marginal_censoring_fit, censoring_models = censoring_models)
}

######################################################################
### get_propensity_scores.R ends here
