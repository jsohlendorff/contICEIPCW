#' @keywords internal
#' @useDynLib contICEIPCW, .registration = TRUE
#' @importFrom survival Surv strata coxph basehaz
#' @importFrom riskRegression predictRisk
#' @importFrom data.table data.table set dcast setkeyv as.data.table copy data.table is.data.table melt rbindlist setnames setorder setorderv setcolorder setkey ":=" ".N" ".SD"
#' @importFrom stats predict glm qnorm pnorm as.formula quasibinomial lm model.matrix dnorm binomial sd gaussian 
NULL

#' @importFrom Rcpp sourceCpp
.onLoad <- function(libname, pkgname) {
  # Ensure Rcpp is initialized (safe even if not strictly required)
  invisible()
}

.datatable.aware <- TRUE

## Deal with idiotic notes from R CMD check about variables in data.table code
utils::globalVariables(c(
  "A", "A_0", "L", "L_0", "Lambda_C_diff",
  "N", "A", "A_0", "L", "L_0",
  "cens_mg", "entrytime", "event", "event_k_prev",
  "event_number", "id", "integrand", "inverse_cumulative_probability_weights",
  "ipw_cum_weight", "last_event", "n_A_events", "n_L_events",
  "new_A", "pseudo_outcome", "protocol_follow", "q_pred_u",
  "q_prediction", "surv", "terminal_time", "time", "time_horizon",
  "time_k_prev", "type", ".", "..baseline_covariates", "q_diff", "q_prediction_prev"
))
