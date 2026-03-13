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
