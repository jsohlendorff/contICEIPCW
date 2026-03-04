### censoring_martingale.R --- 
#----------------------------------------------------------------------
## Author: Johan Sebastian Ohlendorff
## Created: Feb 27 2026 (15:06) 
## Version: 
## Last-Updated: Mar  4 2026 (19:28) 
##           By: Johan Sebastian Ohlendorff
##     Update #: 399
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
censoring_martingale <- function(
  data_censoring,
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
  static_intervention
) {

  ## ------------------------------------------------------------------
  ## 0. Checks
  ## ------------------------------------------------------------------

  if (is.null(grid_size) || grid_size < 1)
    stop("grid_size must be >= 1.")

  if (!marginal_censoring)
    stop("Currently only implemented for marginal censoring.")

  time_k      <- paste0("time_", k)
  time_k_prev <- paste0("time_", k - 1)
  event_k     <- paste0("event_", k)

  ids <- data_at_risk$id

  ## ------------------------------------------------------------------
  ## 1. Check if censoring exists
  ## ------------------------------------------------------------------

  cens_exists <- data_at_risk[
    get(event_k) == "C" & get(time_k) <= time_horizon,
    .N
  ] > 0

  if (!cens_exists) {
    data_at_risk[, cens_mg := 0]
    return(
      merge(
        data_at_risk[, .(id, pseudo_outcome, q_prediction, cens_mg)],
        data[, .(id, inverse_cumulative_probability_weights)],
        by = "id"
      )[, inverse_cumulative_probability_weights :=
          inverse_cumulative_probability_weights *
          (pseudo_outcome - q_prediction + cens_mg)]
    )
  }

  ## ------------------------------------------------------------------
  ## 2. Fit survival model (event hazard)
  ## ------------------------------------------------------------------

  surv_fit <- hazard_fit(
    at_risk_interevent,
    model_hazard,
    outcome_string =
      paste0("Surv(", time_k, ", ", event_k,
             " %in% c(\"A\",\"L\",\"D\",\"Y\"))"),
    covariates = NULL,
    formula_strategy = "additive",
    use_history_of_variables = TRUE,
    lag = lag,
    k = k,
    time_covariates = time_covariates,
    baseline_covariates = baseline_covariates,
    time_variable = time_k
  )

  if (!inherits(surv_fit$fit, "coxph"))
    stop("Only coxph hazard models currently supported.")

  ## ------------------------------------------------------------------
  ## 3. Build time grid
  ## ------------------------------------------------------------------

  if (k == 1)
    data_at_risk[, time_0 := 0]

  tmin <- min(data_at_risk[[time_k_prev]])
  tmax <- min(max(data_at_risk[[time_k]]), time_horizon)
  tmid <- stats::quantile(data_at_risk[[time_k]], probs = 0.05)
    
  if (time_horizon >= tmid) {
      tupper <- seq(tmid, tmax, length.out = grid_size - 1)
      time_grid <- c(tmin, tupper)
  } else {
      time_grid <- seq(tmin, tmax, length.out = grid_size)
  }  

  ## protocol follow indicator
  data_at_risk[, protocol_follow :=
    data[id, inverse_cumulative_probability_weights] > 0]

  ids_follow <- data_at_risk[protocol_follow == TRUE, id]

  ## ------------------------------------------------------------------
  ## 4. Predict q_u over grid
  ## ------------------------------------------------------------------

  preds <- lapply(time_grid, function(u) {
    dt <- data_at_risk[
      get(time_k_prev) < u
    ]

    if (u == min(time_grid)) {
      dt[, q_pred_u := 0]
    } else if (u == time_horizon) {
      dt <- dt[protocol_follow == TRUE]
      dt[, q_pred_u := q_prediction]
    } else {
      dt[, pseudo_outcome_u :=
            pseudo_outcome * (get(time_k) <= u)]

      q_fit <- regression_fit(
        dt,
        model_pseudo_outcome,
        outcome_string = "pseudo_outcome_u",
        use_history_of_variables = TRUE,
        lag = lag,
        k = k,
        time_covariates = time_covariates,
        baseline_covariates = baseline_covariates,
        type = "pseudo_outcome"
      )

       dt[, q_pred_u := q_fit(.SD)]
    }
    dt <- dt[protocol_follow == TRUE]

    dt[, .(
      id,
      time = u,
      q_diff = pmax(q_prediction - q_pred_u, 0)
    )]
  })

  preds <- rbindlist(preds)
  ## NOTE: Need preds_start for linear approximation of q_diff, start point 
  preds_start <- data.table::data.table(id = ids_follow, time = data_at_risk[id %in% ids_follow, get(time_k_prev)], q_diff = data_at_risk[id %in% ids_follow, q_prediction])
  preds <- rbind(preds, preds_start)
  preds[, type := "pred"]
  setkey(preds, id, time)

  ## ------------------------------------------------------------------
  ## 5. Counting process times
  ## ------------------------------------------------------------------
  cens_times <- data_at_risk[
    protocol_follow == TRUE &
      get(event_k) == "C" &
      get(time_k) <= time_horizon,
    .(id, time = get(time_k))
  ][, `:=`(q_diff = NA_real_, type = "counting_process")]

  ## ------------------------------------------------------------------
  ## 6. Censoring cumulative hazard grid (at all terminal event times)
  ## ------------------------------------------------------------------

  censor_times <- data_censoring[
    time <= time_horizon,
    unique(time)
  ]

  cj_dat <- data.table::CJ(time = censor_times, id = ids_follow)
  cj_dat <- merge(
    cj_dat,
    data_at_risk[, .(id,
                     time_k = get(time_k),
                     time_k_prev = get(time_k_prev))],
    by = "id"
  )

  cj_dat <- cj_dat[
    time_k_prev < time & time <= time_k,
    .(id, time)
  ][, `:=`(q_diff = NA_real_, type = "cumhazard")]

  ## ------------------------------------------------------------------
  ## 7. Combine prediction + hazard structure
  ## ------------------------------------------------------------------

  preds <- rbindlist(list(preds, cens_times, cj_dat), fill = TRUE)
  setkey(preds, id, time)

  ## Approximate q_diff with linear interpolation
  preds[, q_diff := zoo::na.approx(q_diff, x = time,
                                   na.rm = FALSE),
        by = id]
  preds <- preds[type != "pred"]

  ## ------------------------------------------------------------------
  ## 8. Compute event survival S
  ## ------------------------------------------------------------------

  preds_surv <- merge(
    preds,
    data_at_risk[, .(id,time_k_prev = get(time_k_prev))],
    by = "id"
  )

  preds_surv[, time := time - time_k_prev] ## Interarrival event form

  preds_surv <- cumulative_hazard_cox(
    surv_fit$fit,
    preds_surv,
    at_risk_interevent
  )

  preds_surv[, surv := exp(-Lambda_minus)]
  preds_surv[, time := time + time_k_prev]
    
  preds <- merge(
    preds,
    preds_surv[, .(id, time, surv, type)],
    by = c("id", "time", "type")
  )

  ## ------------------------------------------------------------------
  ## 9. Compute censoring hazard
  ## ------------------------------------------------------------------

  if (!inherits(marginal_censoring_fit$fit, "coxph"))
    stop("Only coxph marginal censoring supported.")

  preds <- merge(preds, data_at_risk[, c("id", time_k_prev), with = FALSE], by = "id")

  preds <- cumulative_hazard_cox(
    marginal_censoring_fit$fit,
    preds,
    data_censoring[, !"time"],
    time_ref = time_k_prev
  )
  

  setnames(preds,
           c("Lambda", "Lambda_minus"),
           c("Lambda_C", "Lambda_C_minus"))

  preds[, `:=`(
    surv_cens = exp(-Lambda_C_minus),
    Lambda_C_diff = Lambda_C - Lambda_C_minus
  )]

  ## ------------------------------------------------------------------
  ## 10. Build martingale terms
  ## ------------------------------------------------------------------

  preds[, integrand :=
              q_diff / (surv * surv_cens)]

  counting_term <- preds[type == "counting_process",
                         .(mg_counting_term = integrand),
                         by = id]

  lambda_term <- preds[type == "cumhazard",
                       .(mg_lambda_term =
                           sum(Lambda_C_diff * integrand)),
                       by = id]

  mg <- merge(counting_term,
              lambda_term,
              by = "id",
              all = TRUE)
  mg[is.na(mg)] <- 0
  mg[, cens_mg := mg_counting_term - mg_lambda_term]

  ## ------------------------------------------------------------------
  ## 11. Merge back and return IC
  ## ------------------------------------------------------------------  
  out <- merge(
    data_at_risk,
    mg[, .(id, cens_mg)],
    by = "id",
    all.x = TRUE
  )

  out[is.na(cens_mg), cens_mg := 0]

  ic_final <- merge(
    out[, .(id, pseudo_outcome,
            q_prediction, cens_mg)],
    data[, .(id,
             inverse_cumulative_probability_weights)],
    by = "id"
  )

  ic_final[
    , inverse_cumulative_probability_weights :=
        inverse_cumulative_probability_weights *
        (pseudo_outcome - q_prediction + cens_mg)
  ]
}

######################################################################
### censoring_martingale.R ends here
