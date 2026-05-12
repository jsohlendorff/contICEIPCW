### learn_Q.R --- 
#----------------------------------------------------------------------
## Author: Johan Sebastian Ohlendorff
## Created: Mar 13 2026 (18:49) 
## Version: 
## Last-Updated: May 12 2026 (18:13) 
##           By: Johan Sebastian Ohlendorff
##     Update #: 287
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

# Model to use for the outcome regression which returns a prediction function
# for the chosen model.
# Available models are:
# \code{"ranger"},
# \code{"lm"}.
learn_Q <- function(model_type,
                    history_of_variables,
                    data_learn,
                    formula_strategy = "additive",
                    outcome_name = "weight",
                    outcome_string_unweighted = NULL,
                    ipcw_name = NULL,
                    penalize,
                    verbose,
                    k,
                    save_coefficients = FALSE) {
    max_weight <- max(data_learn[[outcome_name]])
    coefficients <- NULL
    if (is.null(max_weight) || is.na(max_weight)) {
        stop("The 'weight' column in data_learn must not be NULL or NA.")
    }
    if (max_weight == 0) {
        predict_fun <- function(data) {
            rep(0, nrow(data))
        }
        return(list(predict_fun = predict_fun, coefficients = "0 values"))
    }
    
    ## String for the formula of the rhs of the formula with the outcome
    if (length(history_of_variables) == 0) {
        history_of_variables_string <- "1"
    } else if (formula_strategy == "additive") {
        history_of_variables_string <- paste(history_of_variables, collapse = "+")
    } else {
        stop("Currently only 'additive' formula strategy is supported.")
    }

    if (verbose) message("Fitting outcome regression model of type: ", model_type, " with formula: ", outcome_name, "_", k, " ~ ", history_of_variables_string)

    if (model_type %in% c("lm", "ipcw_glm_expit", "ipcw_glm_probit")) {
        if (model_type %in% c("ipcw_glm_expit", "ipcw_glm_probit")) {
            ## See Eq. (2) of https://link.springer.com/article/10.1007/s10985-022-09564-6
            ## Can be implemented by fitting a glm with weights;
            ## although we do need censoring survival weights at time $tau$ for that.

            data_learn$out <- data_learn[[outcome_string_unweighted]]
            weights <- data_learn[[ipcw_name]]
            if (model_type == "ipcw_glm_expit") {
                family <- quasibinomial(link = "logit")
            } else if (model_type == "ipcw_glm_probit") {
                family <- binomial(link = "probit")
            }
        } else {
            data_learn$out <- data_learn[[outcome_name]]
            weights <- rep(1, nrow(data_learn))
            family <- stats::gaussian()
        }

        if (!penalize || length(history_of_variables) == 1){ ## do not run penalized regression with one covariate only
              fit <- stats::glm(
                as.formula(paste0("out ~", history_of_variables_string
                )),
                data = data_learn,
                family = family,
                weights = weights
              )
              predict_fun <- function(data) {
                predict(fit, data, type = "response") 
              }
              if (save_coefficients) {
                  coefficients <- coef(fit)
                  coefficients[is.na(coefficients)] <- 0
              }
       } else {
           ## Use Lasso with glmnet
           X <- model.matrix(as.formula(paste0(" ~ ", history_of_variables_string)), data = data_learn)[,-1]
           y <- data_learn[["out"]]
           tryCatch({
               cv_fit <- glmnet::cv.glmnet(X, y, alpha = 1, weights = weights, family = family)
           },
           error = function(e) {
                if (verbose) {
                     warning("glmnet cv.glmnet failed: ", e$message)
                     message("Trying different lambda sequence for glmnet...")
                } 
                lambdas <- glmnet::glmnet(X, y, alpha = 1, weights = weights, family = family)$lambda
                cv_fit <<- glmnet::cv.glmnet(X, y, alpha = 1, weights = weights, lambda = lambdas, family = family)
           })
           
           fit <- glmnet::glmnet(X, y, alpha = 1, lambda = cv_fit$lambda.min, weights = weights, family = family)
              predict_fun <- function(data) {
                X_new <- model.matrix(as.formula(paste0(" ~ ", history_of_variables_string)), data = data)[, -1]
                as.vector(predict(fit, newx = X_new, s = "lambda.min", type = "response"))
              }
           if (save_coefficients) {
               coefs <- coef(fit)
               coefficients <- as.numeric(coefs)
               names(coefficients) <- rownames(coefs)
           }
       }
    } else if (model_type %in% c("oipcw_expit", "oipcw_probit", "nls_expit", "nls_probit")) {
       Y <- data_learn[[outcome_name]]
       X_all <- model.matrix(as.formula(paste0(" ~ ", history_of_variables_string)), data = data_learn)

       ## Remove columns of X which are NA in qr.coef(qr(X), Y);
       ## i.e., columns that would normally be removed in a glm fit
       qr_coef <- qr.coef(qr(X_all),Y)
       X <- X_all[, !is.na(qr_coef), drop = FALSE]

       ## Start with initial parameters corresponding to intercept only model
       intercept <- ifelse(grepl("expit", model_type), logit(mean(Y)), stats::qnorm(mean(Y)))
       beta_init <- c(intercept, rep(0, ncol(X) - 1))

       ## Determine if expit or probit from model_type
       link_function <- ifelse(grepl("expit", model_type), expit, stats::pnorm)

       if (model_type == "oipcw_expit") {
           g <- function(beta, X, Y) {
               eta <- X %*% beta
               as.vector(t(X) %*% (Y - expit(eta)))
           }
       } else if (model_type == "oipcw_probit") {
           g <- function(beta, X, Y) {
               eta <- X %*% beta
               as.vector(t(X) %*% ((Y - stats::pnorm(eta)) * dnorm(eta) / (stats::pnorm(eta) * (1 - stats::pnorm(eta)))))
           }
       } else if (model_type == "nls_expit") {
           g <- function(beta, X, Y) {
               eta <- X %*% beta
               mu <- expit(eta)
               as.vector(t(X) %*% ((Y - mu) * mu * (1 - mu)))
           } 
       } else if (model_type == "nls_probit") {
           g <- function(beta, X, Y) {
               eta <- X %*% beta
               mu <- stats::pnorm(eta)
               as.vector(t(X) %*% ((Y - mu) * dnorm(eta)))
           }
       }
       
       ## Call to fast C++ implementation of the estimating equation solver (ChatGPT)
       ## only use this with OIPCW;
       ## NLS has serious issues with the C++ solver
       tryCatch({
           if (grepl("oipcw", model_type)) {
           fit <- as.vector(suppressWarnings(estimating_equation_cpp(
               X = X,
               Y = Y,
               maxit = 1000,
               tol = 1e-8,
               beta = beta_init,
               offset = rep(0, nrow(X))
           )))
       } else {
           requireNamespace("nleqslv", quietly = TRUE)
           fit <- nleqslv::nleqslv(f = g, x = beta_init, X = X, Y = Y, control = list(maxit = 1000, allowSingular = TRUE))$x
       }},
       error = function(e) {
           warning("The estimating equation solver did not converge: ", e$message)
           fit <<- beta_init
       })

       check_if_failed <- function(fit, g, X, Y, beta_init, verbose,
                                   ignore_large_solution = FALSE, name = "cpp") {
           failed <- FALSE
           warnings_list <- character(0)

           if (is.null(fit) || !is.numeric(fit) || anyNA(fit) || any(!is.finite(fit))) {
               if (verbose)
                   warnings_list <- c(warnings_list,
                                      paste("Solver failed for", name, "(NA / non-finite / invalid fit)."))
               return(list(
                   failed = TRUE,
                   warnings = warnings_list,
                   g_val = NA_real_,
                   fit = beta_init
               ))
           }

           if (any(abs(fit) > 1e2) && !ignore_large_solution) {
               failed <- TRUE
           }

           g_val <- tryCatch(as.numeric(g(fit, X, Y)), error = function(e) NULL)

           if (is.null(g_val) || anyNA(g_val) || any(!is.finite(g_val))) {
               failed <- TRUE
               g_val <- NA_real_
           } else if (any(abs(g_val) > 1e-2)) {
               failed <- TRUE
           }

           list(failed = failed, warnings = warnings_list, g_val = g_val, fit = fit)
       }
       check_fit <- check_if_failed(fit, g, X, Y, beta_init, verbose)
       failed <- check_fit$failed
       warnings_fit <- check_fit$warnings
       fit <- check_fit$fit
       g_val <- check_fit$g_val
       abs_g_val <- ifelse(all(!is.na(g_val)), sum(abs(g_val))/length(g_val), Inf)
       if (verbose) message("C++ solver average absolute value of estimating equation function: ", abs_g_val)
       
       if (grepl("oipcw", model_type) && requireNamespace("nleqslv", quietly = TRUE)) {
           tryCatch({
               requireNamespace("nleqslv", quietly = TRUE)
               fit_nleqslv <- nleqslv::nleqslv(f = g, x = beta_init, X = X, Y = Y, control = list(maxit = 1000, allowSingular = TRUE))$x
           },
           error = function(e) {
               if (verbose) warning("The estimating equation solver did not converge: ", e$message)
               fit_nleqslv <<- beta_init
           })
           check_fit_nleqslv <- check_if_failed(fit_nleqslv, g, X, Y, verbose, ignore_large_solution = TRUE, name = "nleqslv")
           warning_nleqslv <- check_fit_nleqslv$warnings
           failed_nleqslv <- check_fit_nleqslv$failed
           fit_nleqslv <- check_fit_nleqslv$fit
           g_val_nleqslv <- check_fit_nleqslv$g_val
           abs_g_val_nleqslv <- ifelse(all(!is.na(g_val_nleqslv)), sum(abs(g_val_nleqslv))/length(g_val_nleqslv), Inf)
           if (verbose) message("nleqslv solver average absolute value of estimating equation function: ", abs_g_val_nleqslv)
           
           g_val_beta_init <- tryCatch(as.numeric(g(beta_init, X, Y)), error = function(e) NA_real_, warning = function(w) NA_real_)
           abs_g_val_beta_init <- ifelse(all(!is.na(g_val_beta_init)), sum(abs(g_val_beta_init))/length(g_val_beta_init), Inf)
           if (verbose) message("Initial beta average absolute value of estimating equation function: ", abs_g_val_beta_init)
           best_index <- which.min(c(abs_g_val_beta_init, abs_g_val, abs_g_val_nleqslv))
           
           ## if (verbose) message("Winner: ", c("beta_init", "cpp_fit", "nleqslv_fit")[best_index], " with value of estimating equation function: ", round(values[best_index],4))
           fit <- list(beta_init, fit, fit_nleqslv)[[best_index]]
           abs_g_val <- c(abs_g_val_beta_init, abs_g_val, abs_g_val_nleqslv)[best_index]
           if (verbose) message("Best method: ", c("beta_init", "cpp_fit", "nleqslv_fit")[best_index], " with average absolute value of estimating equation function: ", abs_g_val)
           if (abs_g_val> 1e-2) {
               warning("solvers failed to solve the estimating equation. Average absolute value of estimating equation function: ", round(abs_g_val,4))
               if (verbose) {
               message("Possible issues: \n")
               for (w in warnings_fit) {
                   message(w)
               }
               for (w in warning_nleqslv) {
                   message(w)
               }
               message("\n")
           }
           }
       } else if (abs_g_val > 1e-2) {
           warning("solver failed to solve the estimating equation. Average absolute value of estimating equation function: ", round(abs_g_val,4))
       }
       if (verbose){
           message("\n")
       }

       predict_fun <- function(data) {
           X_new <- model.matrix(as.formula(paste0(
               " ~ ", history_of_variables_string
           )), data = data)[, !is.na(qr_coef), drop = FALSE]

           link_function(X_new %*% fit)
       }
       if (save_coefficients) {
           coefficients <- vector("numeric", ncol(X_all))
           coefficients[!is.na(qr_coef)] <- fit
           coefficients[is.na(qr_coef)] <- 0
           names(coefficients) <- colnames(X_all)
       }
   } else {
       ## If flexible, we should pick argmin_(f in cal(F)) sum((Y - g(f(X)))^2, where Y are the outcome weights and g is either expit or probit.
       ## This ensures that predictions will be 0/1-valued.
        formula_w <- paste0(
            outcome_name,
            " ~ ", history_of_variables_string
        )
        predict_fun <- do.call(model_type, list(
                                               character_formula = formula_w,
                                               data = data_learn
                                           ))$predict_fun
   }
   return(list(predict_fun = predict_fun, coefficients = coefficients))
}

######################################################################
### learn_Q.R ends here
