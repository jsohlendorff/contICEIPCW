apply_tmle_update <- function(ic_final,
                              q_prediction,
                              solver = estimating_equation_cpp,
                              verbose = FALSE) {
    has_weight <- ic_final$ipw_cum_weight > 0 & !is.na(ic_final$ipw_cum_weight)
    if (any(has_weight)) {
        Y <- ic_final$pseudo_outcome[has_weight]
        offset <- logit(ic_final$q_prediction[has_weight])
        X <- matrix(1, nrow = sum(has_weight), ncol = 1)
        weights <- ic_final$ipw_cum_weight[has_weight]
        weights <- as.numeric(scale(weights, center = FALSE))
        epsilonhat <- tryCatch(
            solver(
                X = X,
                Y = Y,
                maxit = 1000,
                tol = 1e-8,
                beta = 0,
                offset = offset,
                weights_ = weights
            )[1, 1],
            error = function(e) {
                if (verbose) {
                    message("TMLE update failed: ", e$message)
                }
                0
            }
        )
    } else {
        epsilonhat <- 0
    }

    if (is.null(epsilonhat) || length(epsilonhat) != 1 || !is.finite(epsilonhat)) {
        warning("TMLE update failed to compute estimating equation value; No TMLE update performed for this iteration.")
        q_new <- ic_final$q_prediction
    } else {
        q_new <- expit(epsilonhat + logit(ic_final$q_prediction))
    }

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

    list(ic_final = ic_final, q_prediction = q_prediction, epsilonhat = epsilonhat)
}
