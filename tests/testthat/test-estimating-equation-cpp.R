test_that("estimating_equation_cpp validates input dimensions", {
    X <- matrix(1, nrow = 5, ncol = 1)
    Y <- rep(0.5, 5)
    offset <- rep(0, 5)

    expect_error(
        contICEIPCW:::estimating_equation_cpp(
            X = X,
            Y = Y[-1],
            beta = 0,
            offset = offset
        ),
        "Length of Y must match number of rows in X"
    )

    expect_error(
        contICEIPCW:::estimating_equation_cpp(
            X = X,
            Y = Y,
            beta = 0,
            offset = offset[-1]
        ),
        "Length of offset must match number of rows in X"
    )

    expect_error(
        contICEIPCW:::estimating_equation_cpp(
            X = X,
            Y = Y,
            beta = 0,
            offset = offset,
            weights_ = rep(1, 4)
        ),
        "Length of weights_ must match number of rows in X"
    )
})

test_that("estimating_equation_cpp accepts vector weights", {
    X <- matrix(1, nrow = 5, ncol = 1)
    Y <- seq(0.1, 0.5, length.out = 5)
    offset <- rep(0, 5)
    weights <- as.numeric(scale(seq_len(5), center = FALSE))

    fit <- contICEIPCW:::estimating_equation_cpp(
        X = X,
        Y = Y,
        beta = 0,
        offset = offset,
        weights_ = weights
    )

    expect_type(fit, "double")
    expect_equal(length(fit), 1)
    expect_true(is.finite(fit))
})

test_that("estimating_equation_r matches estimating_equation_cpp", {
    set.seed(17)
    X <- cbind(
        1,
        scale(rnorm(80)),
        scale(rnorm(80))
    )
    beta <- c(-0.2, 0.4, -0.3)
    offset <- rnorm(nrow(X), sd = 0.15)
    Y <- contICEIPCW:::expit(as.vector(X %*% beta) + offset)

    fit_cpp <- contICEIPCW:::estimating_equation_cpp(
        X = X,
        Y = Y,
        beta = rep(0, ncol(X)),
        offset = offset,
        maxit = 1000,
        tol = 1e-10
    )
    fit_r <- contICEIPCW:::estimating_equation_r(
        X = X,
        Y = Y,
        beta = rep(0, ncol(X)),
        offset = offset,
        maxit = 1000,
        tol = 1e-10
    )

    expect_equal(dim(fit_r), dim(fit_cpp))
    expect_equal(as.vector(fit_r), as.vector(fit_cpp), tolerance = 1e-8)
})

test_that("estimating_equation_r matches estimating_equation_cpp with weights", {
    set.seed(18)
    X <- cbind(
        1,
        scale(rnorm(60)),
        scale(rnorm(60))
    )
    beta <- c(0.15, -0.25, 0.35)
    offset <- rnorm(nrow(X), sd = 0.1)
    Y <- contICEIPCW:::expit(as.vector(X %*% beta) + offset)
    weights <- seq(0.25, 2.5, length.out = nrow(X))

    fit_cpp <- contICEIPCW:::estimating_equation_cpp(
        X = X,
        Y = Y,
        beta = rep(0, ncol(X)),
        offset = offset,
        weights_ = weights,
        maxit = 1000,
        tol = 1e-10
    )
    fit_r <- contICEIPCW:::estimating_equation_r(
        X = X,
        Y = Y,
        beta = rep(0, ncol(X)),
        offset = offset,
        weights_ = weights,
        maxit = 1000,
        tol = 1e-10
    )

    expect_equal(dim(fit_r), dim(fit_cpp))
    expect_equal(as.vector(fit_r), as.vector(fit_cpp), tolerance = 1e-8)
})

test_that("TMLE update gives the same result with R and C++ solvers", {
    ic_final <- data.table::data.table(
        id = 1:8,
        ipw_cum_weight = c(1, 2, 0, NA, 3, 1.5, 0.5, 4),
        pseudo_outcome = c(0, 0.1, 0.2, 0.3, 0.6, 0.7, 0.8, 1),
        q_prediction = c(0.08, 0.12, 0.25, 0.35, 0.55, 0.65, 0.72, 0.9)
    )
    q_prediction <- ic_final[, .(id, q_prediction_prev = q_prediction)]

    tmle_cpp <- contICEIPCW:::apply_tmle_update(
        ic_final = data.table::copy(ic_final),
        q_prediction = data.table::copy(q_prediction),
        solver = contICEIPCW:::estimating_equation_cpp
    )
    tmle_r <- contICEIPCW:::apply_tmle_update(
        ic_final = data.table::copy(ic_final),
        q_prediction = data.table::copy(q_prediction),
        solver = contICEIPCW:::estimating_equation_r
    )

    expect_equal(tmle_r$epsilonhat, tmle_cpp$epsilonhat, tolerance = 1e-8)
    expect_equal(
        tmle_r$ic_final$q_prediction,
        tmle_cpp$ic_final$q_prediction,
        tolerance = 1e-8
    )
    expect_equal(
        tmle_r$q_prediction$q_prediction_prev,
        tmle_cpp$q_prediction$q_prediction_prev,
        tolerance = 1e-8
    )
})
