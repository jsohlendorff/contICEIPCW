apply_tmle_update <- function(ic_final,
                              q_prediction,
                              verbose = FALSE) {
    has_weight <- ic_final$ipw_cum_weight > 0 & !is.na(ic_final$ipw_cum_weight)
    if (any(has_weight)) {
        pseudo_outcome <- ic_final$pseudo_outcome[has_weight]
        offset <- logit(ic_final$q_prediction[has_weight])
        X <- matrix(1, nrow = sum(has_weight), ncol = 1)
        ipw <- ic_final$ipw_cum_weight[has_weight]
        ipw <- as.numeric(scale(ipw, center = FALSE))
        epsilonhat <- tryCatch(
            estimating_equation_cpp(
                X = X,
                Y = pseudo_outcome,
                maxit = 1000,
                tol = 1e-8,
                beta = 0,
                offset = offset,
                weights_ = ipw
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
    estimating_equation_tmle_step <- function(epsilon, ipw, pseudo_outcome,q_pred) {
        as.vector(t(ipw) %*% (pseudo_outcome - expit(logit(q_pred) + epsilon)))
    }
    estimating_equation_tmle_step_alternate <- function(epsilon, ipw, pseudo_outcome,q_pred) {
        as.vector(t(ipw) %*% (pseudo_outcome - expit(logit(q_pred) + epsilon * ipw)))
    }

    check_if_failed <- function(epsilon, estimating_equation_fn, ipw, pseudo_outcome, q_pred, tol = 1e-3) {
        valid <- TRUE
        failed <- FALSE
        if (is.null(epsilonhat) || length(epsilonhat) != 1 || !is.finite(epsilonhat)) {
            valid <- FALSE
            failed <- TRUE
        } 

        ## Check if the estimating equation value is solved
        g_val <- estimating_equation_fn(epsilon, ipw, pseudo_outcome, q_pred)
        if (abs(g_val) > tol) {
            failed <- TRUE
        }
        return(list(failed = failed, g_val = g_val, valid = valid))
    }
    check_result <- check_if_failed(epsilonhat, estimating_equation_tmle_step, ic_final$ipw_cum_weight, ic_final$pseudo_outcome, ic_final$q_prediction, tol = 1e-5)
    failed <- check_result$failed
    g_val <- check_result$g_val
    valid_cpp <- check_result$valid

    ## Retry with nleqslv if the initial solution did not solve the estimating equation
    if (failed) {
        requireNamespace("nleqslv")
        ## Try solving both and selecting the best 
        tryCatch(
        {
            ## ipw <- ic_final$ipw_cum_weight[has_weight]
            q_pred <- ic_final$q_prediction[has_weight]
            epsilonhat_nleqslv <- nleqslv::nleqslv(
                                               x = 0,
                                               fn = function(epsilon) estimating_equation_tmle_step(epsilon, ipw, pseudo_outcome, q_pred),
                                               control = list(maxit = 1000, allowSingular = TRUE)
                                           )$x
            epsilonhat_nleqslv_alternate <- nleqslv::nleqslv(
                                                         x = 0,
                                                         fn = function(epsilon) estimating_equation_tmle_step_alternate(epsilon, ipw, pseudo_outcome, q_pred),
                                                         control = list(maxit = 1000, allowSingular = TRUE)
                                                     )$x
            ## Optimizer
            epsilonhat_nleqslv_minimizer <- stats::optimize(
                lower = -0.5, upper = 0.5,
                f = function(epsilon) abs(estimating_equation_tmle_step(epsilon, ipw, pseudo_outcome, q_pred))
            )$minimum
            
            check_nleqslv <- check_if_failed(epsilonhat_nleqslv, estimating_equation_tmle_step, ipw, pseudo_outcome, q_pred, tol = 1e-5)
            check_nleqslv_alternate <- check_if_failed(epsilonhat_nleqslv_alternate, estimating_equation_tmle_step_alternate, ipw, pseudo_outcome, q_pred, tol = 1e-5)
            check_nleqslv_minimizer <- check_if_failed(epsilonhat_nleqslv_minimizer, estimating_equation_tmle_step, ipw, pseudo_outcome, q_pred, tol = 1e-5)
            check_zero <- check_if_failed(0, estimating_equation_tmle_step, ipw, pseudo_outcome, q_pred, tol = 1e-5)

            ## Pick best solution among the three candidates (amonh those with valid estimating equation values)
            candidates <- data.frame(
                epsilon = c(epsilonhat, epsilonhat_nleqslv, epsilonhat_nleqslv_alternate, epsilonhat_nleqslv_minimizer, 0),
                g_val = c(g_val, check_nleqslv$g_val, check_nleqslv_alternate$g_val, check_nleqslv_minimizer$g_val, check_zero$g_val),
                valid = c(valid_cpp, check_nleqslv$valid, check_nleqslv_alternate$valid, check_nleqslv_minimizer$valid, check_zero$valid)
            )
            valid_candidates <- candidates[candidates$valid, ]
            if (nrow(valid_candidates) == 0) {
                warning("All TMLE update attempts failed to solve the estimating equation; No TMLE update performed for this iteration.")
                epsilonhat <- 0
            } else {
                best_candidate <- valid_candidates[which.min(abs(valid_candidates$g_val)), ]
                epsilonhat <- best_candidate$epsilon
                if (verbose) {
                    message("TMLE update: Selected epsilon = ", epsilonhat, " with estimating equation value = ", best_candidate$g_val)
                }
                if (best_candidate$g_val > 1e-3) {
                    warning("TMLE update resulted in a solution with estimating equation value = ", best_candidate$g_val, ". Solution may be degenerate. ")
                }
            }
        },
        error = function(e) {
            warning("TMLE update with nleqslv failed: ", e$message, " No TMLE update performed for this iteration.")
            epsilonhat <<- 0
        })
    }

    ## Define q_new based on the selected epsilonhat
    ## If the best candidate is alt use expit(epsilonhat * ic_final$ipw_cum_weight + logit(ic_final$q_prediction)) instead of expit(epsilonhat + logit(ic_final$q_prediction))
    if (exists("best_candidate") && !is.null(best_candidate) && best_candidate$epsilon == epsilonhat_nleqslv_alternate) {
        q_new <- expit(epsilonhat * ic_final$ipw_cum_weight + logit(ic_final$q_prediction))
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
