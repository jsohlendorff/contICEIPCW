### censoring_martingale.R --- 
#----------------------------------------------------------------------
## Author: Johan Sebastian Ohlendorff
## Created: Feb 27 2026 (15:06) 
## Version: 
## Last-Updated: Mar  3 2026 (09:50) 
##           By: Johan Sebastian Ohlendorff
##     Update #: 197
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
censoring_martingale <- function(data_censoring,
                                 at_risk_before_time_horizon,
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
                                 static_intervention) {
    event_k<-time_k<-time_0<-protocol_follow<-inverse_cumulative_probability_weights<-id <- time_k_minus <- pseudo_outcome_f <- pseudo_outcome <- q_pred_u <- q_prediction <- q_diff <- type <- time <- init <- Lambda <- Lambda_minus <- surv <- Lambda_0 <- Lambda_C <- Lambda_C_diff <- surv_cens <- integrand <- mg_lambda_term <- cens_mg <- mg_counting_term <- . <- NULL
    if (!is.null(grid_size) && grid_size > 0) {
        if (!marginal_censoring) {
            stop("multiple_ice is only supported for marginal censoring hazard (for now).")
        }
        ## If grid size is specified, we need to create a grid of time points
        ## and interpolate the predictions on that grid.

        utimes <- unique(sort(data_censoring$time))
        min_utimes <- min(diff(c(0, utimes))) / 2
        exists_censored_event_k <- nrow(at_risk_before_time_horizon[event_k == "C" & time_k <= time_horizon, env = list(
                                                                                                                 event_k = paste0("event_", k),
                                                                                                                 time_k = paste0("time_", k)
                                                                                                             )]) > 0
        ## Learn survival function
        if (exists_censored_event_k) {
            learn_surv <- hazard_fit(at_risk_interevent,
                                     model_hazard,
                                     outcome_string = paste0("Surv(time_", k, ", event_", k, " %in% c(\"A\", \"L\", \"D\", \"Y\"))"),
                                     covariates = NULL,
                                     formula_strategy = "additive",
                                     use_history_of_variables = TRUE,
                                     lag = lag,
                                     k = k,
                                     time_covariates = time_covariates,
                                     baseline_covariates = baseline_covariates)
            
            utimes_surv <- unique(sort(at_risk_interevent[[paste0("time_", k)]]))
            min_utimes_surv <- min(diff(c(0, utimes_surv))) / 2

            if (k == 1) {
                at_risk_before_time_horizon[, time_0 := 0]
            }
            time_grid_min <- stats::quantile(at_risk_before_time_horizon[[paste0("time_", k-1)]], probs = 0)
            time_grid_mid <- stats::quantile(at_risk_before_time_horizon[[paste0("time_", k)]], probs = c(0.05,1))
            ## Check if time_horizon < time_grid_mid[1]; if not do the following
            if (time_horizon >= time_grid_mid[1]) {
                time_grid_upper <- seq(time_grid_mid[1], min(time_grid_mid[2], time_horizon), length.out = grid_size - 1)
                time_grid <- c(time_grid_min, time_grid_upper)
            } else {
                time_grid_max <- stats::quantile(at_risk_before_time_horizon[[paste0("time_", k)]], probs = 1)
                ## Uniform grid between 0 and time_horizon
                time_grid <- seq(time_grid_min, pmin(time_grid_max, time_horizon), length.out = grid_size)
            }                        
            
            ids <- at_risk_before_time_horizon$id
            at_risk_before_time_horizon[, protocol_follow := data[ids, inverse_cumulative_probability_weights] > 0]
            ids_follow <- at_risk_before_time_horizon[protocol_follow == TRUE, id]
            preds <- list()
            for (v in seq_len(grid_size)) {
                ##print(paste0("Processing time ", time_grid[v]))
                ## Should probably do it with at_risk_before_u
                if (v > 1) {
                    at_risk_before_u <- at_risk_before_time_horizon[time_k_minus < time_grid[v], env = list(time_k_minus = paste0("time_", k - 1))]
                    at_risk_before_u[, pseudo_outcome_f := pseudo_outcome * 1 * (time_k <= time_grid[v]), env = list(time_k = paste0("time_", k))]
                    if (time_grid[v] == time_horizon) {
                        at_risk_before_u <- at_risk_before_u[protocol_follow == TRUE]
                        at_risk_before_u[, q_pred_u := q_prediction]
                    } else {
                        q_hat_u <- regression_fit(
                            at_risk_before_u,
                            model_pseudo_outcome,
                            outcome_string = "pseudo_outcome_f",
                            use_history_of_variables = TRUE,
                            lag = lag,
                            k = k,
                            time_covariates = time_covariates,
                            baseline_covariates = baseline_covariates,
                            type = "pseudo_outcome"
                        )
                    }
                    at_risk_before_u <- at_risk_before_u[protocol_follow == TRUE]
                    at_risk_before_u[, q_pred_u := predict_intervention(.SD, k-1, q_hat_u, static_intervention)] 
                } else {
                    at_risk_before_u <- at_risk_before_time_horizon[time_k_minus < time_grid[v], env = list(time_k_minus = paste0("time_", k - 1))]
                    at_risk_before_u <- at_risk_before_u[protocol_follow == TRUE]
                    at_risk_before_u[, q_pred_u := 0]
                }
                data_u <- at_risk_before_u[, c("id", "q_pred_u", "q_prediction"), with = FALSE]
                ## connect data_u with data_at_risk_time_horizon ids not available; set their q_pred_u to 0
                ids_time_horizon <- setdiff(ids, at_risk_before_u$id)
                if (length(ids_time_horizon) > 0) {
                    ## if (length(ids_time_horizon) != nrow(data_u)) browser()
                    ## if (k == 2) 
                    data_u <- rbind(
                        data_u,
                        data.frame(id = ids_time_horizon, q_pred_u = 0, q_prediction = at_risk_before_time_horizon[id %in% ids_time_horizon, q_prediction])
                    )
                    
                }
                data_u$time <- time_grid[v]
                data_u$q_diff <-  data_u$q_prediction - data_u$q_pred_u
                preds[[v]] <- data_u[, c("id", "time", "q_diff"), with = FALSE]
            }
            preds <- rbindlist(preds)
            preds[q_diff < 0, q_diff := 0]
            preds[, type := "pred"]
            ## Counting process term
            ## Determine which observations we need for N^c and S^c, S, q_time_horizon q_u
            dat_cens_times <- at_risk_before_time_horizon[protocol_follow & time_k <= time_horizon & event_k == "C", env = list(time_k = paste0("time_", k),
                                                                                                                                event_k = paste0("event_", k))][, c("id", paste0("time_", k)), with = FALSE]
            dat_cens_times[, c("q_diff", "type") := list(NA, "counting_process")]
            colnames(dat_cens_times)[2] <- "time"
            
            ## Compensator term
            ## Determine which observations we need for N^c and S^c, S, q_time_horizon q_u
            data_censoring_times_k <- data_censoring[time <= time_horizon, time]
            min_cens_prev <- min(at_risk_before_time_horizon[[paste0("time_", k - 1)]])
            data_censoring_times_k <- data_censoring_times_k[data_censoring_times_k >= min_cens_prev]
            cj_dat <- data.table::CJ(time = data_censoring_times_k, id = ids_follow)
            cj_dat <- merge(cj_dat, at_risk_before_time_horizon[, c("id", paste0("time_", k), paste0("time_", k-1)), with = FALSE], by = "id")
            cj_dat <- cj_dat[time_k_minus < time & time <= time_k & time <= time_horizon,
                             env = list(time_k = paste0("time_", k), time_k_minus = paste0("time_", k - 1))]
            cj_dat <- cj_dat[, c("id", "time")]
            ## Message which ids are ids but not in cj_dat
            ids_in_cj_dat <- unique(cj_dat$id)
            ids_not_in_cj_dat <- setdiff(ids, ids_in_cj_dat)
            ## if (length(ids_not_in_cj_dat) > 0) {
            ##     message("The following ids are not in cj_dat: ", paste(ids_not_in_cj_dat, collapse = ", "))
            ## }
            cj_dat[, q_diff := NA]
            cj_dat[, type := "cumhazard"]

            ## Combine preds, dat_cens_times, and cj_dat
            preds <- rbindlist(list(preds, dat_cens_times, cj_dat))
            setkeyv(preds, c("id", "time"))

            ## Approximate places with NAs using interpolation, i.e., the "approx" function in R, by id
            preds[, q_diff := zoo::na.approx(q_diff, x = time, na.rm = FALSE), by = id]
            ## Remove superflous rows (again?)
            preds <- merge(preds, at_risk_before_time_horizon[, c("id", paste0("time_", k), paste0("time_", k-1)), with = FALSE], by = "id")
            preds <- preds[time_k_minus <= time & time <= time_k & time <= time_horizon, env = list(
                                                                                             time_k = paste0("time_", k),
                                                                                             time_k_minus = paste0("time_", k - 1)
                                                                                         )]
            preds <- preds[type != "pred", c("id", "time", "q_diff", "type")]
            ## Add to preds the time of the k-1'th event by id
            ## This is needed for diff Cumhazard at the next time point.
            preds_temp <- at_risk_before_time_horizon[, c("id", paste0("time_", k-1)), with = FALSE][, c("q_diff", "type") := list(NA, "cumhazard")]
            colnames(preds_temp)[2] <- "time"
            preds_temp[, init := TRUE]
            preds[, init := FALSE]
            preds <- rbind(preds, preds_temp)
            setkeyv(preds, c("id", "time"))

            if (!inherits(learn_surv$fit, "coxph")){
                stop("Not implemented!")
            } else {
                ## Because pred's may not be on interevent form; we need to make sure that it is
                preds_mod <- copy(preds)
                ## Add to preds the time of the k-1'th event by id
                preds_mod <- merge(preds_mod, at_risk_before_time_horizon[, c("id", paste0("time_", k-1)), with = FALSE], by = "id")
                ## Subtract time_k-1 from time in preds_mod to get the time since the last event
                preds_mod[, time_org := time]
                preds_mod[, time := time - get(paste0("time_", k-1))]
                setnames(preds_mod, paste0("time_", k-1), "time_k_minus")
                preds_mod <- merge(preds_mod, at_risk_interevent, by = "id")
                preds_mod <- cumulative_hazard_cox(learn_surv$fit, preds_mod)
                ## cumulative_hazard_cox_variant(learn_surv, preds_mod, at_risk_interevent[id %in% preds$id])
                ## Change NAs to zero
                ## preds_mod[is.na(Lambda), Lambda := 0]
                ## preds_mod[is.na(Lambda_minus), Lambda_minus := 0]
                preds_mod[, surv := exp(-Lambda_minus)]
                ## preds_mod[, c("Lambda", "Lambda_minus") := list(NULL, NULL)]
                ## preds_mod[, time := time + get(paste0("time_", k-1))]
                ## preds_mod[, time_k_minus := NULL, env = list(time_k_minus = paste0("time_", k - 1))]
                preds_mod <- preds_mod[, c("id", "time_org", "surv", "type")]
                setnames(preds_mod, "time_org", "time")
                preds <- unique(merge(preds, preds_mod[, c("id", "time","surv", "type")], by = c("id", "time", "type")))
                ## Because pred's may not be on interevent form; we need to make sure that it is 
                ## preds_mod <- preds
                ## ## Add to preds the time of the k-1'th event by id
                ## preds_mod <- merge(preds_mod, at_risk_before_time_horizon[, c("id", paste0("time_", k-1)), with = FALSE], by = "id")
                ## ## Subtract time_k-1 from time in preds_mod to get the time since the last event
                ## preds_mod[, time := time - get(paste0("time_", k-1))]
                ## preds_mod <- cumulative_hazard_cox_variant(learn_surv, preds_mod, at_risk_interevent[id %in% preds$id])
                ## ## Change NAs to zero
                ## preds_mod[is.na(Lambda), Lambda := 0]
                ## preds_mod[is.na(Lambda_minus), Lambda_minus := 0]
                ## preds_mod[, surv := exp(-Lambda_minus)]
                ## preds_mod[, c("Lambda", "Lambda_minus") := list(NULL, NULL)]
                ## preds_mod[, time := time + get(paste0("time_", k-1))]
                ## preds_mod[, time_k_minus := NULL, env = list(time_k_minus = paste0("time_", k - 1))]
                ## preds <- preds_mod
            }
            if (!inherits(marginal_censoring_fit$fit, "coxph")){
                ## preds <- pred_and_merge(preds, marginal_censoring_fit$fit, data_censoring[id %in% preds$id], unique(preds$time), unique(preds$time), "surv_cens")
                stop("*Not implemented!")
            } else {
                ## Find previous event time for each id in preds; i.e., time_{k-1}
                preds <- merge(preds, at_risk_before_time_horizon[, c("id", paste0("time_", k-1)), with = FALSE], by = "id")
                setnames(preds, paste0("time_", k-1), "time_k_minus")
                data_censoring_temp <- data_censoring[id %in% preds$id]
                data_censoring_temp[, time := NULL]
                ## also merge with the data_censoring to get the covariates for the censoring model
                preds <- merge(preds, data_censoring_temp)
                preds_new <- cumulative_hazard_cox(marginal_censoring_fit$fit, preds, time_ref = "time_k_minus")
                setnames(preds_new, c("Lambda", "Lambda_minus"), c("Lambda_C", "Lambda_C_minus"))
                preds <- merge(preds, preds_new, by = c("id", "time"))
                preds[, surv_cens := exp(-Lambda_C_minus)]
                preds[, Lambda_C_diff := Lambda_C - Lambda_C_minus]
            }
            preds <- preds[init == FALSE]
            ## preds[, Lambda_C := -log(surv_cens)]
            ## preds[, Lambda_C_diff := Lambda_C - data.table::shift(Lambda_C, fill = 0), by = id]
            preds[, integrand := q_diff * 1 / surv * 1 / surv_cens]
            ## Counting process term.
            preds_counting <- preds[type == "counting_process"]
            preds_counting <- preds_counting[, c("id", "integrand"), with = FALSE]
            counting_terms <- data.table(id = setdiff(ids, preds_counting$id), integrand = 0)
            preds_counting <- rbind(preds_counting, counting_terms)
            setnames(preds_counting, "integrand", "mg_counting_term")

            ## Lambda term.
            preds_lambda <- preds[type == "cumhazard"][, .(mg_lambda_term = sum(Lambda_C_diff * integrand)), by = id]
            ## Combine the counting process term and the Lambda term
            preds_mg <- merge(preds_counting, preds_lambda, by = "id", all.x = TRUE)
            ## If any NA
            if (any(is.na(preds_mg$mg_lambda_term))) {
                ids_na <- preds_mg[is.na(mg_lambda_term), id]
                ## Warn if ids_na and ids_not_in_cj_dat are not the same
                if (length(ids_not_in_cj_dat) > 0 && length(ids_na) > 0 && !all(ids_na %in% ids_not_in_cj_dat)) {
                    message("The following ids have NA values in mg_lambda_term: ", paste(ids_na, collapse = ", "))
                    message("The following ids are not in cj_dat: ", paste(ids_not_in_cj_dat, collapse = ", "))
                    warning("NA values in MGc, setting them to zero.")
                }
                preds_mg[is.na(mg_lambda_term), mg_lambda_term := 0]
            }
            preds_mg[, cens_mg := mg_counting_term - mg_lambda_term]
            ## Merge with at_risk_before_time_horizon
            at_risk_before_time_horizon <- merge(at_risk_before_time_horizon, preds_mg[, c("id", "cens_mg")], by = "id")
        } else {
            at_risk_before_time_horizon[, cens_mg := 0]
        }
        ic_final <- merge(at_risk_before_time_horizon[, c("pseudo_outcome", "q_prediction", "id", "cens_mg")], data[, c("inverse_cumulative_probability_weights", "id")], by = "id")
        return(ic_final[, inverse_cumulative_probability_weights := inverse_cumulative_probability_weights * (pseudo_outcome - q_prediction + cens_mg)])
    } else {
        stop("Grid size must be specified for conservative=FALSE.")
    }
}

cumulative_hazard_cox <- function(fit, data, time_ref = NULL) {
    hazard<-hazard_minus<-Lambda<-exp_lp <- Lambda_minus <- NULL
    ## Find exp(LP); i.e., exponential of linear predictor
    exp_lp_dt <- data.table(id = data$id)
    exp_lp_dt$exp_lp <- predict(fit,
                                newdata = data,
                                type = "risk",
                                reference = "zero")
    exp_lp_dt <- unique(exp_lp_dt)
    ## Baseline cumulative hazard Lambda_0^x (T_j) for all j
    base_hazard <- as.data.table(survival::basehaz(fit, centered = FALSE))
    base_hazard[, hazard_prev := c(0, hazard[-.N])]

    if (!is.null(time_ref)) {
        min_time <- min(data[, time_ref, with = FALSE])
        data <- data.table::melt(data, id.vars = "id", measure.vars = c("time", time_ref), value.name = "time")
        data <- unique(data)
        setkeyv(data, c("id", "time"))
    } else {
        min_time <- 0
    }
    ## Add 0 row to base_hzazrd at min_time
    if (min_time < base_hazard$time[1]) {
        base_hazard <- rbindlist(list(data.table(hazard = 0, time = min_time, hazard_prev = 0), base_hazard))
    }
    base_hazard_matched_exact <- base_hazard[data, on = "time", nomatch = 0]
    base_hazard_matched_exact[, exact_match := TRUE]

    base_hazard_all_matched <- base_hazard[data, on = "time", roll = TRUE]
    base_hazard_all_matched[, exact_match := FALSE]
    
    ## Merge/roll forward and calculate cumulative hazard function
    data <- base_hazard_all_matched[base_hazard_matched_exact, exact_match := TRUE, on = "time"]
    data[exact_match == FALSE, hazard_minus := hazard]
    data[exact_match == TRUE, hazard_minus := hazard_prev]
    data[, exact_match := NULL]
    data[, hazard_prev := NULL]
    
    if (!is.null(time_ref)) {
        data_time_ref <- data[variable == time_ref]
        data_time_ref[time == 0, "hazard" := 0]
        data_time_ref[, c("hazard_minus", "variable", "time"):= NULL]
        setnames(data_time_ref, "hazard", "hazard_time_ref")
        data <- data[variable != time_ref]
        data[, variable := NULL]
        data <- merge(data, data_time_ref, by = "id", all.x = TRUE)
    } else {
        data[, hazard_time_ref := 0]
    }
    data <- merge(data, exp_lp_dt, by = "id")
    data[, Lambda := exp_lp * (hazard - hazard_time_ref)]
    data[, Lambda_minus := exp_lp * (hazard_minus - hazard_time_ref)]
    data[, c("exp_lp", "hazard", "hazard_minus", "hazard_time_ref") := NULL]
    data
}


######################################################################
### censoring_martingale.R ends here
