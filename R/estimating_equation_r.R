estimating_equation_r <- function(X,
                                  Y,
                                  beta,
                                  offset,
                                  weights_ = NULL,
                                  maxit = 100L,
                                  tol = 1e-8,
                                  verbose = FALSE) {
    X <- as.matrix(X)
    Y <- as.numeric(Y)
    beta <- as.numeric(beta)
    offset <- as.numeric(offset)

    p <- ncol(X)
    n <- nrow(X)

    if (length(beta) != p) {
        stop("Length of beta must match number of columns in X")
    }
    if (length(Y) != n) {
        stop("Length of Y must match number of rows in X")
    }
    if (length(offset) != n) {
        stop("Length of offset must match number of rows in X")
    }

    has_weights <- !is.null(weights_)
    if (has_weights) {
        weights <- as.numeric(weights_)
        if (length(weights) != n) {
            stop("Length of weights_ must match number of rows in X")
        }
    }

    for (iter in seq_len(maxit)) {
        eta <- as.vector(X %*% beta) + offset
        mu <- expit(eta)
        w <- mu * (1 - mu)

        if (has_weights) {
            score <- as.vector(crossprod(X, weights * (Y - mu)))
            jacobian <- -crossprod(X, X * (weights * w))
        } else {
            score <- as.vector(crossprod(X, Y - mu))
            jacobian <- -crossprod(X, X * w)
        }

        step <- tryCatch(
            solve(jacobian, score),
            error = function(e) NULL
        )
        if (is.null(step)) {
            warning("Linear system could not be solved, stopping iteration")
            return(matrix(unname(beta), ncol = 1))
        }
        step <- as.numeric(step)

        if (verbose) {
            message("step is ", paste(step, collapse = " "))
        }

        step_factor <- 1
        beta_new <- beta - step
        while (!all(is.finite(beta_new)) && step_factor > 1e-6) {
            step_factor <- step_factor * 0.5
            beta_new <- beta - step_factor * step
            if (verbose) {
                message("beta_new is ", paste(beta_new, collapse = " "))
            }
        }

        if (max(abs(beta_new - beta)) < tol) {
            return(matrix(unname(beta_new), ncol = 1))
        }

        beta <- beta_new
    }

    warning("Did not converge")
    matrix(unname(beta), ncol = 1)
}
