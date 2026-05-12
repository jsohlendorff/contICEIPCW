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

    score_fn <- function(beta) {
        eta <- as.vector(X %*% beta) + offset
        mu <- expit(eta)
        if (has_weights) {
            as.vector(crossprod(X, weights * (Y - mu)))
        } else {
            as.vector(crossprod(X, Y - mu))
        }
    }

    jacobian_fn <- function(beta) {
        eta <- as.vector(X %*% beta) + offset
        mu <- expit(eta)
        w <- mu * (1 - mu)
        if (has_weights) {
            -crossprod(X, X * (weights * w))
        } else {
            -crossprod(X, X * w)
        }
    }

    fit <- tryCatch(
        nleqslv::nleqslv(
            x = beta,
            fn = score_fn,
            jac = jacobian_fn,
            method = "Newton",
            control = list(
                maxit = maxit,
                xtol = tol,
                ftol = tol,
                trace = if (verbose) 1L else 0L
            )
        ),
        error = function(e) e
    )

    if (inherits(fit, "error")) {
        warning("Root solver failed: ", conditionMessage(fit))
        return(matrix(unname(beta), ncol = 1))
    }

    converged <- identical(fit$termcd, 1L) ||
        isTRUE(max(abs(fit$fvec)) < tol)
    if (!converged) {
        warning(if (!is.null(fit$message)) fit$message else "Did not converge")
    }

    matrix(unname(as.numeric(fit$x)), ncol = 1)
}
