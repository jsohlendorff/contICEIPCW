#' @keywords internal
#' @useDynLib contICEIPCW, .registration = TRUE
#' @importFrom ggplot2 aes aes_string element_blank element_line element_rect geom_errorbar geom_line geom_point geom_ribbon ggplot labs guide_legend guides scale_colour_manual scale_color_continuous scale_fill_manual scale_linetype_manual scale_y_continuous theme theme_bw "%+replace%"  unit xlab  ylab
#' @importFrom survival Surv strata coxph basehaz
#' @importFrom prodlim Hist prodlim
#' @importFrom riskRegression predictRisk
#' @importFrom data.table data.table set dcast setkeyv as.data.table copy data.table is.data.table melt rbindlist setnames setorder setorderv setcolorder setkey ":=" ".N" ".SD"
#' @importFrom stats predict glm qnorm pnorm as.formula quasibinomial lm model.matrix dnorm binomial sd 
NULL

#' @importFrom Rcpp sourceCpp
.onLoad <- function(libname, pkgname) {
  # Ensure Rcpp is initialized (safe even if not strictly required)
  invisible()
}

.datatable.aware = TRUE
