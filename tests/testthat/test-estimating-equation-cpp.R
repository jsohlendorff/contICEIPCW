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
