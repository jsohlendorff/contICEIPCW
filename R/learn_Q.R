### learn_Q.R --- 
#----------------------------------------------------------------------
## Author: Johan Sebastian Ohlendorff
## Created: Mar 13 2026 (18:49) 
## Version: 
## Last-Updated: Apr  1 2026 (11:50) 
##           By: Johan Sebastian Ohlendorff
##     Update #: 108
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
# \code{"quasibinomial"},
# \code{"scaled_quasibinomial"},
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
                    k) {
    max_weight <- max(data_learn[[outcome_name]])
    if (is.null(max_weight) || is.na(max_weight)) {
        stop("The 'weight' column in data_learn must not be NULL or NA.")
    }
    if (max_weight == 0) {
        predict_fun <- function(data) {
            rep(0, nrow(data))
        }
        return(predict_fun)
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

    if (model_type %in% c("quasibinomial", "scaled_quasibinomial", "lm", "ipcw_glm_expit", "ipcw_glm_probit")) {
        if (grepl("quasibinomial", model_type)) {
            if (model_type == "quasibinomial") {
                scale <- 1
            } else if (model_type == "scaled_quasibinomial") {
                scale <- max_weight
                penalize <- FALSE
            }
            set(
                data_learn,
                j = "out",
                value = data_learn[[outcome_name]] / scale
            )
            weights <- rep(1, nrow(data_learn))
            family <- quasibinomial()
        } else if (model_type %in% c("ipcw_glm_expit", "ipcw_glm_probit")) {
            ## See Eq. (2) of https://link.springer.com/article/10.1007/s10985-022-09564-6
            ## Can be implemented by fitting a glm with weights;
            ## although we do need censoring survival weights at time $tau$ for that.

            scale <- 1
            data_learn$out <- data_learn[[outcome_string_unweighted]]
            weights <- data_learn[[ipcw_name]]
            if (model_type == "ipcw_glm_expit") {
                family <- quasibinomial(link = "logit")
            } else if (model_type == "ipcw_glm_probit") {
                family <- binomial(link = "probit")
            }
        } else {
            scale <- 1
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
                predict(fit, data, type = "response") * scale
              }
       } else {
           ## Use Lasso with glmnet
           X <- model.matrix(as.formula(paste0(" ~ ", history_of_variables_string)), data = data_learn)
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
                X_new <- model.matrix(as.formula(paste0(" ~ ", history_of_variables_string)), data = data)
                as.vector(predict(fit, newx = X_new, s = "lambda.min", type = "response"))
              }
       }
    } else if (model_type %in% c("oipcw_expit", "oipcw_probit", "nls_expit", "nls_probit")) {
       Y <- data_learn[[outcome_name]]
       X <- model.matrix(as.formula(paste0(" ~ ", history_of_variables_string)), data = data_learn)

       ## Remove columns of X which are NA in qr.coef(qr(X), Y);
       ## i.e., columns that would normally be removed in a glm fit
       qr_coef <- qr.coef(qr(X),Y)
       X <- X[, !is.na(qr_coef), drop = FALSE]

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
               model_type = model_type,
               maxit = 1000,
               tol = 1e-8,
               beta = beta_init,
               solve_opts = "force_approx",
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


       check_if_failed <- function(fit, g, X, Y, verbose, ignore_large_solution = FALSE, name = "cpp") {
           failed <- FALSE
           warnings_list <- character(0)

           ## Check for NAs or NULL in solution
           if (any(is.na(fit)) || any(is.null(fit))) {
               if (verbose) warnings_list <- c(warnings_list, paste("The estimating equation solver did not converge for method", name, "and returned NA or NULL values."))
               fit <- beta_init
               failed <- TRUE
           }
           
           ## Check if the solution is very large or if the estimating equation does not seem to be solved
           if (any(abs(fit) > 1e2) && !ignore_large_solution) {
               if (verbose) {
                   warnings_list <- c(warnings_list,
                                      paste0("The solution of the estimating equation solver is very large for ", name ,": ", paste(round(fit, 4), collapse = ", "), "\n",
                                             "This may indicate non-convergence or a problem with the estimating equation.")) 
               }
               failed <- TRUE
           }

           g_val <- g(fit, X, Y)
           if (any(abs(g_val) > 1e-2)){
               if (verbose) {
                   warnings_list <- c(warnings_list,
                                      paste0("The estimating equation does not seem to be solved which may indicate non-convergence for ", name, ". \n",
                                             "The value of the estimating equation at the solution is: ", paste(round(g_val, 4), collapse = ", ")))
               }
               failed <- TRUE
           }
           return(list(failed = failed, warnings = warnings_list, g_val = g_val))
       }
       check_fit <- check_if_failed(fit, g, X, Y, verbose)
       failed <- check_fit$failed
       warnings_fit <- check_fit$warnings
       
       if (failed){
           if (grepl("oipcw", model_type) && requireNamespace("nleqslv", quietly = TRUE)) {
               if (verbose) message("Trying to solve the estimating equation with nleqslv...")
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
               
               if (!failed_nleqslv) {
                   fit <- fit_nleqslv
                   if (verbose) message("nleqslv successfully solved the estimating equation.")
               } else {
                   warning("cpp and nleqslv solvers failed to solve the estimating equation. Picking the best solution among, beta_init, fit, and fit_nleqslv based on the value of the estimating equation.")
                   values <- sapply(list(beta_init, fit, fit_nleqslv), function(beta) sum(g(beta, X, Y)^2))
                   best_index <- which.min(values)
                   ## best_index_format
                   
                   if (verbose) message("Winner: ", c("beta_init", "cpp_fit", "nleqslv_fit")[best_index], " with value of estimating equation function: ", round(values[best_index],4))
                   fit <- list(beta_init, fit, fit_nleqslv)[[best_index]]
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
           } else {
               warning("The estimating equation may not be solved. \n ")
               for (w in warnings_fit) {
                   message(w)
               }
               message("\n")
           }
       }

       ## Checks and if nleqslv is avail
       if (failed && grepl("oipcw", model_type)) {
           if (requireNamespace("nleqslv", quietly = TRUE)) {
               
           }            
       }

       predict_fun <- function(data) {
           X_new <- model.matrix(as.formula(paste0(
               " ~ ", history_of_variables_string
           )), data = data)[, !is.na(qr_coef), drop = FALSE]

           link_function(X_new %*% fit)
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
    predict_fun
}

######################################################################
### learn_Q.R ends here
