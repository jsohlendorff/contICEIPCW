## Utility functions for the package

## Simulate from an exponential proportional hazards model
rexponential_proportional_hazard <- function(n, rate, eta) {
    u <- stats::runif(n)
    (-log(u) / (rate * exp(eta)))
}

expit <- function(x) 1 / (1 + exp(-x))
logit <- function(p) log(p / (1 - p))
