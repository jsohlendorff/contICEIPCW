### prepare_data.R --- 
#----------------------------------------------------------------------
## Author: Johan Sebastian Ohlendorff
## Created: Mar  4 2026 (19:33) 
## Version: 
## Last-Updated: Mar  6 2026 (14:33) 
##           By: Johan Sebastian Ohlendorff
##     Update #: 36
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

#' @title Prepare data for continuous time analysis
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
#' @param max_time_horizon A numeric value representing the maximum time horizon considered for the analysis.
#' @param time_covariates A character vector of column names in \code{data} that are
#'   treated as time-varying covariates. Must include values of time-varying covariates at baseline.
#' @param baseline_covariates A character vector of column names in \code{data} that are
#'   considered baseline (time-invariant) covariates. Must include treatment and time-varying covariates.
#' @param last_non_terminal_event Optional numeric indicating the last nonterminal event number to consider
#'   in the outcome.
#' @param marginal_censoring Logical; if \code{TRUE}, assumes censoring depends only on baseline covariates.
#' @export
#'
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

prepare_data <- function(data,
                         max_time_horizon,
                         time_covariates,
                         baseline_covariates,
                         marginal_censoring = TRUE,
                         last_non_terminal_event = NULL) {
    event_number <- id <- ic <- pseudo_outcome <- survival_censoring_k <- event_k <- time_k <- inverse_cumulative_probability_weights <- inverse_cumulative_probability_weights_k_prev <- ipw <- ipw_k <- pred_0 <- estimate <- g_formula_estimate <- . <- NULL
    ## Check user input
    check_input(baseline_covariates, time_covariates, data, max_time_horizon)

    ## Get timevarying data and baseline data and add event number by id
    timevarying_data <- data$timevarying_data[, event_number := seq_len(.N), by = id]
    baseline_data <- data$baseline_data

    ## If last event number not provided,
    ## select last event number adaptively because the iterative
    ## regressions may not have sufficient data to fit the models for later events.
    ## NOTE: Modifies data.
    select_last_event_out <- select_last_event(timevarying_data, max_time_horizon, last_non_terminal_event)
    timevarying_data <- select_last_event_out$timevarying_data
    last_event <- select_last_event_out$last_event
    
    censoring_info_result <- censoring_info(timevarying_data,
                                            baseline_data,
                                            max_time_horizon,
                                            marginal_censoring)
    is_censored <- censoring_info_result$is_censored
    data_marginal_censoring <- censoring_info_result$data_marginal_censoring
    
    ## Convert the data from long format to wide format
    wide_data <- widen_continuous_data(timevarying_data,
                                       baseline_data,
                                       time_covariates)

    out <- list(
        wide_data = wide_data,
        is_censored = is_censored,
        data_marginal_censoring = data_marginal_censoring,
        last_event = last_event,
        time_covariates = time_covariates,
        baseline_covariates = baseline_covariates,
        marginal_censoring = marginal_censoring
    )
    class(out) <- "prepare_data_continuous"
    return(out)
}

######################################################################
### prepare_data.R ends here
