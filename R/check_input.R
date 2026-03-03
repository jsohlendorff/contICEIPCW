check_input <- function(baseline_covariates,
                        time_covariates,
                        data,
                        time_horizon) {
    id<-time<-N<-.<-..<-NULL
    ## TODO: Need to more thorougly check user input.
    ## Check the baseline covariates and time covariates are not empty and character strings
    if (length(baseline_covariates) == 0 || !is.character(baseline_covariates)) {
        stop("baseline_covariates must be a non-empty character vector.")
    }

    if (length(time_covariates) == 0 || !is.character(time_covariates)) {
        stop("time_covariates must be a non-empty character vector.")
    }

    ## Check that data is a list with two data frames, called timevarying_data
    ## and baseline_data
    if (!is.list(data) || !all(c("timevarying_data", "baseline_data") %in% names(data))) {
        stop("Data must be a list containing two data frames: 'timevarying_data' and 'baseline_data'.")
    }

    ## Check that timevarying_data and baseline_data are data.tables
    if (!data.table::is.data.table(data$timevarying_data) || !data.table::is.data.table(data$baseline_data)) {
        stop("Both 'timevarying_data' and 'baseline_data' must be data.tables.")
    }

    ## Check that id, A_0, baseline_covariates are in baseline_data
    if (!all(c("id", "A_0", baseline_covariates) %in% names(data$baseline_data))) {
        stop("Baseline data must contain 'id', 'A_0', and all baseline covariates.")
    }

    ## Check that id, time, event, A, time_covariates are in timevarying_data
    if (!all(c("id", "time", "event", "A", time_covariates) %in% names(data$timevarying_data))) {
        stop("Time-varying data must contain 'id', 'time', 'event', 'A', and all time covariates.")
    }

    ## Check that time_horizon is a positive numeric value
    if (!is.numeric(time_horizon) || length(time_horizon) != 1 || time_horizon <= 0) {
        stop("tau must be a positive numeric value.")
    }

    ## Check that event is a factor with that includes A, L, C, Y, D, and end_of_study
    if (!is.factor(data$timevarying_data$event) || !all(levels(data$timevarying_data$event) %in% c("A", "L", "C", "Y", "D", "tauend"))) {
        stop("The 'event' column in time-varying data must be a factor with levels 'A', 'L', 'C', 'Y', 'D', and 'tauend'.")
    }

    ## Check that time is numeric and non-negative
    if (!is.numeric(data$timevarying_data$time) || any(data$timevarying_data$time < 0)) {
        stop("The 'time' column in time-varying data must be numeric and non-negative.")
    }

    ## Check that the variables specified in time_covariates and baseline_covariates do not contain NULLs or NAs
    if (any(sapply(data$timevarying_data[, .SD, .SDcols = time_covariates], function(x) any(is.null(x) | is.na(x))))) {
        stop("Time-varying covariates must not contain NULL or NA values.")
    }
    if (any(sapply(data$baseline_data[, .SD, .SDcols = baseline_covariates], function(x) any(is.null(x) | is.na(x))))) {
        stop("Baseline covariates must not contain NULL or NA values.")
    }

    ## Check for ties in event times by id
    ties_check <- data$timevarying_data[, .N, by = .(id, time)][N > 1]
    if (nrow(ties_check) > 0) {
        stop("There are ties in event times for some ids. Please ensure that each id has unique event times.")
    }
}

######################################################################
### check_input.R ends here
