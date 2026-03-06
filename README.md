
<!-- README.md is generated from README.Rmd. Please edit that file -->

# contICEIPCW

<!-- badges: start -->

[![R-CMD-check](https://github.com/jsohlendorff/contICEIPCW/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jsohlendorff/contICEIPCW/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The goal of contICEIPCW is to provide an implementation of the ICE-IPCW
algorithm for longitudinal causal inference in continuous-time with
targeted learning for time-to-event outcomes.

## Installation

You can install the development version of contICEIPCW from
[GitHub](https://github.com/) with:

``` r
devtools::install_github("jsohlendorff/contICEIPCW")
```

## Example

``` r
library(contICEIPCW)
#> Loading required package: data.table
set.seed(15)
data_continuous <- simulate_continuous_time_data(
  n = 1000,
  uncensored = FALSE,
  no_competing_events = FALSE,
  baseline_rate_list = list(
    A = 0.005,
    L = 0.001,
    C = 0.0008,
    Y = 0.0001,
    D = 0.00015
  )
)
prep_data <- prepare_data(
  data = data_continuous,
  max_time_horizon = 720,
  time_covariates = c("A", "L"),
  baseline_covariates = c("age", "A_0", "L_0"),
  marginal_censoring = TRUE
)
propensity_score_data <- propensity_scores(
  prepared_data = prep_data,
  model_treatment = "learn_glm_logistic",
  model_hazard = "learn_coxph",
  verbose = TRUE)
result <- debias_ice_ipcw(
  prepared_data = propensity_score_data,
  time_horizon = 720,
  model_pseudo_outcome = "lm",
  model_hazard = "learn_coxph",
  conservative = TRUE,
  verbose = TRUE
)
```
