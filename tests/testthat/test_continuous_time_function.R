library(testthat)

test_that("test continuous time function (uncensored)", {
    library(data.table)

    set.seed(34)
    
    # Simulate continuous time data with continuous and irregular event times
    data_continuous <- simulate_continuous_time_data(
        n = 1000,
        no_competing_events = TRUE,
        uncensored = TRUE
    )

    prep_data <- prepare_data(
        data = data_continuous,
        max_time_horizon = 720,
        time_covariates = c("A", "L"),
        baseline_covariates = c("age", "A_0", "L_0")
    )
    altered_data <- propensity_scores(
        prepared_data = prep_data,
        model_treatment = "learn_glm_logistic",
        model_hazard = NULL
    )

    # Run debiased ICE-IPCW procedure
    result <- debias_ice_ipcw(
        prepared_data = altered_data,
        time_horizon = 720,
        model_pseudo_outcome = "quasibinomial",
        model_hazard = NULL,
        conservative = TRUE,
        verbose = FALSE
    )

    correct_result <- data.table::data.table(
                                      estimate = c(0.282948697909507),
                                      se = c(0.0165150186058971),
                                      lower = c(0.250579261441949),
                                      upper = c(0.315318134377065),
                                      ice_ipcw_estimate = c(0.283170820624823),
                                      ipw = c(0.282975855447292)
                                  )

    expect_true(all.equal(result, correct_result, tolerance = 1e-8))
})

test_that("test continuous time function (censored; conservative)", {
    library(survival)
    library(data.table)
    library(prodlim)
    library(riskRegression)

    set.seed(34)
    # Simulate continuous time data with continuous and irregular event times
    data_continuous <- simulate_continuous_time_data(
        n = 1000,
        no_competing_events = TRUE,
        uncensored = FALSE
    )

    # Run debiased ICE-IPCW procedure
    prep_data <- prepare_data(
        data = data_continuous,
        max_time_horizon = 720,
        time_covariates = c("A", "L"),
        baseline_covariates = c("age", "A_0", "L_0"),
        marginal_censoring = FALSE
    )
    altered_data <- propensity_scores(
        prepared_data = prep_data,
        model_treatment = "learn_glm_logistic",
        model_hazard = "learn_coxph"
    )

    # Run debiased ICE-IPCW procedure
    result <- debias_ice_ipcw(
        prepared_data = altered_data,
        time_horizon = 720,
        model_pseudo_outcome = "scaled_quasibinomial",
        model_hazard = "learn_coxph",
        conservative = TRUE,
        verbose = FALSE
    )

    correct_result <- data.table::data.table(
                                      estimate = c(0.270426500711251),
                                      se = c(0.0167728967294177),
                                      lower = c(0.237551623121592),
                                      upper = c(0.30330137830091),
                                      ice_ipcw_estimate = c(0.271534334306832),
                                      ipw = c(0.269324346522728)
                                  )

    expect_true(all.equal(result, correct_result, tolerance = 1e-8))
})

test_that("test continuous time function (censored; conservative; marginal_censoring_hazard)", {
    library(survival)
    library(data.table)
    library(prodlim)
    library(riskRegression)

    set.seed(34)
    # Simulate continuous time data with continuous and irregular event times
    data_continuous <- simulate_continuous_time_data(
        n = 1000,
        no_competing_events = TRUE,
        uncensored = FALSE
    )

    # Run debiased ICE-IPCW procedure
    prep_data <- prepare_data(
        data = data_continuous,
        max_time_horizon = 720,
        time_covariates = c("A", "L"),
        baseline_covariates = c("age", "A_0", "L_0"),
        marginal_censoring = TRUE
    )
    altered_data <- propensity_scores(
        prepared_data = prep_data,
        model_treatment = "learn_glm_logistic",
        model_hazard = "learn_coxph"
    )

    # Run debiased ICE-IPCW procedure
    result <- debias_ice_ipcw(
        prepared_data = altered_data,
        time_horizon = 720,
        model_pseudo_outcome = "scaled_quasibinomial",
        model_hazard = "learn_coxph",
        conservative = TRUE,
        verbose = FALSE
    )

    correct_result <- data.table::data.table(
                                      estimate = 0.2704025237816331,
                                      se = 0.016772235506564332,
                                      lower = 0.23752894218876702,
                                      upper = 0.30327610537449917,
                                      ice_ipcw_estimate = 0.2714325697499465,
                                      ipw = 0.2693019050719549
                                  )
    expect_true(all.equal(result, correct_result, tolerance = 1e-8))
})

test_that("test continuous time function (censored; non_conservative; multiple ice)", {
    library(survival)
    library(data.table)
    library(prodlim)
    library(riskRegression)

    set.seed(34)
    # Simulate continuous time data with continuous and irregular event times
    data_continuous <- simulate_continuous_time_data(
        n = 1000,
        no_competing_events = TRUE,
        uncensored = FALSE
    )

    prep_data <- prepare_data(
        data = data_continuous,
        max_time_horizon = 720,
        time_covariates = c("A", "L"),
        baseline_covariates = c("age", "A_0", "L_0"),
        marginal_censoring = TRUE
    )
    altered_data <- propensity_scores(
        prepared_data = prep_data,
        model_treatment = "learn_glm_logistic",
        model_hazard = "learn_coxph"
    )

    # Run debiased ICE-IPCW procedure
    result <- debias_ice_ipcw(
        prepared_data = altered_data,
        time_horizon = 720,
        model_pseudo_outcome = "scaled_quasibinomial",
        model_hazard = "learn_coxph",
        conservative = FALSE,
        grid_size = 10,
        verbose = FALSE
    )

    # dpasta(result)
    correct_result <- data.table::data.table(
                                      estimate = 0.2701755590505636,
                                      se = 0.016741446624264904,
                                      lower = 0.2373623236670044,
                                      upper = 0.30298879443412285,
                                      ice_ipcw_estimate = 0.2714325697499465,
                                      ipw = 0.2693019050719549
                                  )

    expect_true(all.equal(result, correct_result, tolerance = 1e-8))
})

test_that("error when time-varying covariates contain NAs", {
    library(data.table)

    set.seed(34)
    # Simulate continuous time data with continuous and irregular event times
    data_continuous <- simulate_continuous_time_data(
        n = 1000,
        no_competing_events = TRUE,
        uncensored = TRUE
    )
    data_continuous$timevarying_data[event == "tauend", L := NA]
    expect_error(
        prepare_data(
        data = data_continuous,
        max_time_horizon = 720,
        time_covariates = c("A", "L"),
        baseline_covariates = c("age", "A_0", "L_0"),
        marginal_censoring = TRUE
    ),
        "Time-varying covariates must not contain NULL or NA values."
    )
})

test_that("error when time-varying covariates contain ties", {
    library(data.table)

    set.seed(34)
    # Simulate continuous time data with continuous and irregular event times
    data_continuous <- simulate_continuous_time_data(
        n = 1000,
        no_competing_events = TRUE,
        uncensored = TRUE
    )
    data_continuous$timevarying_data[id == "2", time := 5]
    expect_error(
        prepare_data(
        data = data_continuous,
        max_time_horizon = 720,
        time_covariates = c("A", "L"),
        baseline_covariates = c("age", "A_0", "L_0"),
        marginal_censoring = TRUE
    ),
        "There are ties in event times for some ids. Please ensure that each id has unique event times"
    )
})

test_that("semiTMLE option", {
    library(data.table)

    set.seed(34)
    # Simulate continuous time data with continuous and irregular event times
    data_continuous <- simulate_continuous_time_data(
        n = 1000,
        no_competing_events = TRUE,
        uncensored = TRUE
    )
        # Simulate continuous time data with continuous and irregular event times
    data_continuous <- simulate_continuous_time_data(
        n = 1000,
        no_competing_events = TRUE,
        uncensored = TRUE
    )

    prep_data <- prepare_data(
        data = data_continuous,
        max_time_horizon = 720,
        time_covariates = c("A", "L"),
        baseline_covariates = c("age", "A_0", "L_0")
    )
    altered_data <- propensity_scores(
        prepared_data = prep_data,
        model_treatment = "learn_glm_logistic",
        model_hazard = NULL
    )

    # Run debiased ICE-IPCW procedure
    expect_no_error(result <- debias_ice_ipcw(
        prepared_data = altered_data,
        time_horizon = 720,
        model_pseudo_outcome = "quasibinomial",
        model_hazard = NULL,
        conservative = TRUE,
        verbose = FALSE,
        semi_tmle = TRUE
    ))
})

test_that("test continuous time function (uncensored; competing risks)", {
    library(data.table)

    set.seed(34)
    # Simulate continuous time data with continuous and irregular event times
    data_continuous <- simulate_continuous_time_data(
        n = 1000,
        no_competing_events = FALSE,
        uncensored = TRUE
    )

    # Run debiased ICE-IPCW procedure
        prep_data <- prepare_data(
        data = data_continuous,
        max_time_horizon = 720,
        time_covariates = c("A", "L"),
        baseline_covariates = c("age", "A_0", "L_0")
    )
    altered_data <- propensity_scores(
        prepared_data = prep_data,
        model_treatment = "learn_glm_logistic",
        model_hazard = NULL
    )

    # Run debiased ICE-IPCW procedure
    result <- debias_ice_ipcw(
        prepared_data = altered_data,
        time_horizon = 720,
        model_pseudo_outcome = "quasibinomial",
        model_hazard = NULL,
        conservative = TRUE,
        verbose = FALSE
    )

    # library(datapasta)
    # dpasta(result)
    correct_result <- data.table::data.table(
                                      estimate = 0.2596497654234689,
                                      se = 0.016206822172967386,
                                      lower = 0.2278843939644528,
                                      upper = 0.291415136882485,
                                      ice_ipcw_estimate = 0.2604858194566298,
                                      ipw = 0.2593272643831233
                                  )

    expect_true(all.equal(result, correct_result, tolerance = 1e-8))
})

test_that("test continuous time function (censored; conservative; competing risks)", {
    library(survival)
    library(data.table)
    library(prodlim)
    library(riskRegression)

    set.seed(34)
    # Simulate continuous time data with continuous and irregular event times
    data_continuous <- simulate_continuous_time_data(
        n = 1000,
        no_competing_events = FALSE,
        uncensored = FALSE
    )

    # Run debiased ICE-IPCW procedure
    prep_data <- prepare_data(
        data = data_continuous,
        max_time_horizon = 720,
        time_covariates = c("A", "L"),
        baseline_covariates = c("age", "A_0", "L_0"),
        marginal_censoring = FALSE
    )
    altered_data <- propensity_scores(
        prepared_data = prep_data,
        model_treatment = "learn_glm_logistic",
        model_hazard = "learn_coxph"
    )

    # Run debiased ICE-IPCW procedure
    result <- debias_ice_ipcw(
        prepared_data = altered_data,
        time_horizon = 720,
        model_pseudo_outcome = "scaled_quasibinomial",
        model_hazard = "learn_coxph",
        conservative = TRUE,
        verbose = FALSE
    )
    
    correct_result <- data.table::data.table(
                                      estimate = 0.28509305096789817,
                                      se = 0.01687112579034439,
                                      lower = 0.25202564441882314,
                                      upper = 0.3181604575169732,
                                      ice_ipcw_estimate = 0.28498047354582223,
                                      ipw = 0.2852953695578807
                                  )

    expect_true(all.equal(result, correct_result, tolerance = 1e-8))
})

test_that("test continuous time function (censored; competing events; conservative; marginal_censoring)", {
    library(survival)
    library(data.table)
    library(prodlim)
    library(riskRegression)

    set.seed(34)
    # Simulate continuous time data with continuous and irregular event times
    data_continuous <- simulate_continuous_time_data(
        n = 1000,
        no_competing_events = FALSE,
        uncensored = FALSE
    )

    prep_data <- prepare_data(
        data = data_continuous,
        max_time_horizon = 720,
        time_covariates = c("A", "L"),
        baseline_covariates = c("age", "A_0", "L_0"),
        marginal_censoring = TRUE
    )
    altered_data <- propensity_scores(
        prepared_data = prep_data,
        model_treatment = "learn_glm_logistic",
        model_hazard = "learn_coxph"
    )

    # Run debiased ICE-IPCW procedure
    result <- debias_ice_ipcw(
        prepared_data = altered_data,
        time_horizon = 720,
        model_pseudo_outcome = "scaled_quasibinomial",
        model_hazard = "learn_coxph",
        conservative = TRUE,
        verbose = FALSE
    )

    correct_result <- data.table::data.table(
                                      estimate = 0.28476294065658975,
                                      se = 0.01684744267744183,
                                      lower = 0.25174195300880375,
                                      upper = 0.31778392830437574,
                                      ice_ipcw_estimate = 0.2846643753705404,
                                      ipw = 0.28495738884128674
                                  )
    expect_true(all.equal(result, correct_result, tolerance = 1e-8))
})

test_that("test continuous time function (censored; competing events; conservative; marginal_censoring, multiple_ice", {
    library(survival)
    library(data.table)
    library(prodlim)
    library(riskRegression)

    set.seed(34)
    # Simulate continuous time data with continuous and irregular event times
    data_continuous <- simulate_continuous_time_data(
        n = 1000,
        no_competing_events = FALSE,
        uncensored = FALSE
    )

    prep_data <- prepare_data(
        data = data_continuous,
        max_time_horizon = 720,
        time_covariates = c("A", "L"),
        baseline_covariates = c("age", "A_0", "L_0"),
        marginal_censoring = TRUE
    )
    altered_data <- propensity_scores(
        prepared_data = prep_data,
        model_treatment = "learn_glm_logistic",
        model_hazard = "learn_coxph"
    )

    # Run debiased ICE-IPCW procedure
    result <- debias_ice_ipcw(
        prepared_data = altered_data,
        time_horizon = 720,
        model_pseudo_outcome = "scaled_quasibinomial",
        model_hazard = "learn_coxph",
        conservative = FALSE,
        grid_size = 15,
        verbose = FALSE
    )

    correct_result <- data.table::data.table(
                                      estimate = 0.28419787544646624,
                                      se = 0.01681271790106795,
                                      lower = 0.25124494836037303,
                                      upper = 0.31715080253255945,
                                      ice_ipcw_estimate = 0.2846643753705404,
                                      ipw = 0.28495738884128674
                                  )
    expect_true(all.equal(result, correct_result, tolerance = 1e-8))    
})

test_that("test continuous time function (censored; competing events; oicpw_expit)", {
    library(survival)
    library(data.table)
    library(prodlim)
    library(riskRegression)

    set.seed(34)
    # Simulate continuous time data with continuous and irregular event times
    data_continuous <- simulate_continuous_time_data(
        n = 1000,
        no_competing_events = FALSE,
        uncensored = FALSE
    )

    # Run debiased ICE-IPCW procedure
    
    prep_data <- prepare_data(
        data = data_continuous,
        max_time_horizon = 720,
        time_covariates = c("A", "L"),
        baseline_covariates = c("age", "A_0", "L_0"),
        marginal_censoring = TRUE
    )
    altered_data <- propensity_scores(
        prepared_data = prep_data,
        model_treatment = "learn_glm_logistic",
        model_hazard = "learn_coxph"
    )

    # Run debiased ICE-IPCW procedure
    result <- debias_ice_ipcw(
        prepared_data = altered_data,
        time_horizon = 720,
        model_pseudo_outcome = "oipcw_expit",
        model_hazard = "learn_coxph",
        conservative = TRUE,
        verbose = FALSE
    )

    correct_result <- data.table::data.table(
                                      estimate = 0.28476394511718484,
                                      se = 0.01684741059880183,
                                      lower = 0.25174302034353324,
                                      upper = 0.31778486989083643,
                                      ice_ipcw_estimate = 0.28466287527342254,
                                      ipw = 0.28495738884128674
                                  )
    expect_true(all.equal(result, correct_result, tolerance = 1e-8))
})

test_that("test continuous time function (censored; competing events; oicpw_probit)", {
    library(survival)
    library(data.table)
    library(prodlim)
    library(riskRegression)
    
    set.seed(34)
    # Simulate continuous time data with continuous and irregular event times
    data_continuous <- simulate_continuous_time_data(
        n = 1000,
        no_competing_events = FALSE,
        uncensored = FALSE
    )

    # Run debiased ICE-IPCW procedure
    
    prep_data <- prepare_data(
        data = data_continuous,
        max_time_horizon = 720,
        time_covariates = c("A", "L"),
        baseline_covariates = c("age", "A_0", "L_0"),
        marginal_censoring = TRUE
    )
    altered_data <- propensity_scores(
        prepared_data = prep_data,
        model_treatment = "learn_glm_logistic",
        model_hazard = "learn_coxph"
    )

    # Run debiased ICE-IPCW procedure
    result <- debias_ice_ipcw(
        prepared_data = altered_data,
        time_horizon = 720,
        model_pseudo_outcome = "oipcw_probit",
        model_hazard = "learn_coxph",
        conservative = TRUE,
        verbose = FALSE
    )

    correct_result <- data.table::data.table(
                                      estimate = 0.2847569509671835,
                                      se = 0.016846434902731817,
                                      lower = 0.25173793855782917,
                                      upper = 0.31777596337653785,
                                      ice_ipcw_estimate = 0.284527778014968,
                                      ipw = 0.28495738884128674
                                  )
    expect_true(all.equal(result, correct_result, tolerance = 1e-8))
})

test_that("test continuous time function (censored; competing events; nls_expit", {
    library(survival)
    library(data.table)
    library(prodlim)
    library(riskRegression)
    
    set.seed(34)
    # Simulate continuous time data with continuous and irregular event times
    data_continuous <- simulate_continuous_time_data(
        n = 1000,
        no_competing_events = FALSE,
        uncensored = FALSE
    )

    
    prep_data <- prepare_data(
        data = data_continuous,
        max_time_horizon = 720,
        time_covariates = c("A", "L"),
        baseline_covariates = c("age", "A_0", "L_0"),
        marginal_censoring = TRUE
    )
    altered_data <- propensity_scores(
        prepared_data = prep_data,
        model_treatment = "learn_glm_logistic",
        model_hazard = "learn_coxph"
    )

    # Run debiased ICE-IPCW procedure
    result <- debias_ice_ipcw(
        prepared_data = altered_data,
        time_horizon = 720,
        model_pseudo_outcome = "nls_expit",
        model_hazard = "learn_coxph",
        conservative = TRUE,
        verbose = FALSE
    )

    correct_result <- data.table::data.table(
                                      estimate = 0.284825969619289,
                                      se = 0.016845604816601346,
                                      lower = 0.2518085841787504,
                                      upper = 0.31784335505982764,
                                      ice_ipcw_estimate = 0.2857275251295884,
                                      ipw = 0.28495738884128674
                                  )
    expect_true(all.equal(result, correct_result, tolerance = 1e-8))
})

test_that("test continuous time function (censored; competing events; nls_probit)", {
    library(survival)
    library(data.table)
    library(prodlim)
    library(riskRegression)

    set.seed(34)
    # Simulate continuous time data with continuous and irregular event times
    data_continuous <- simulate_continuous_time_data(
        n = 1000,
        no_competing_events = FALSE,
        uncensored = FALSE
    )

    
    prep_data <- prepare_data(
        data = data_continuous,
        max_time_horizon = 720,
        time_covariates = c("A", "L"),
        baseline_covariates = c("age", "A_0", "L_0"),
        marginal_censoring = TRUE
    )
    altered_data <- propensity_scores(
        prepared_data = prep_data,
        model_treatment = "learn_glm_logistic",
        model_hazard = "learn_coxph"
    )

    # Run debiased ICE-IPCW procedure
    result <- debias_ice_ipcw(
        prepared_data = altered_data,
        time_horizon = 720,
        model_pseudo_outcome = "nls_probit",
        model_hazard = "learn_coxph",
        conservative = TRUE,
        verbose = FALSE
    )

    correct_result <- data.table::data.table(
                                      estimate = 0.2848191584885638,
                                      se = 0.0168447482523257,
                                      lower = 0.25180345191400544,
                                      upper = 0.3178348650631222,
                                      ice_ipcw_estimate = 0.2855452132800225,
                                      ipw = 0.28495738884128674
                                  )
    expect_true(all.equal(result, correct_result, tolerance = 1e-8))
})

test_that("test continuous time function (censored; competing events; ipcw_glm_expit)", {
    library(survival)
    library(data.table)
    library(prodlim)
    library(riskRegression)

    set.seed(34)
    # Simulate continuous time data with continuous and irregular event times
    data_continuous <- simulate_continuous_time_data(
        n = 1000,
        no_competing_events = FALSE,
        uncensored = FALSE
    )

    
    prep_data <- prepare_data(
        data = data_continuous,
        max_time_horizon = 720,
        time_covariates = c("A", "L"),
        baseline_covariates = c("age", "A_0", "L_0"),
        marginal_censoring = TRUE
    )
    altered_data <- propensity_scores(
        prepared_data = prep_data,
        model_treatment = "learn_glm_logistic",
        model_hazard = "learn_coxph"
    )

    # Run debiased ICE-IPCW procedure
    result <- debias_ice_ipcw(
        prepared_data = altered_data,
        time_horizon = 720,
        model_pseudo_outcome = "ipcw_glm_expit",
        model_hazard = "learn_coxph",
        conservative = TRUE,
        verbose = FALSE
    )

    correct_result <- data.table::data.table(
  estimate = 0.28477536762863775,
  se = 0.0168471119505276,
  lower = 0.25175502820560364,
  upper = 0.31779570705167187,
  ice_ipcw_estimate = 0.28426232748420016,
  ipw = 0.28495738884128674
)
    expect_true(all.equal(result, correct_result, tolerance = 1e-8))
})

test_that("test continuous time function (censored; competing events; lm, penalize)", {
    library(survival)
    library(data.table)
    library(prodlim)
    library(riskRegression)

    set.seed(34)
    # Simulate continuous time data with continuous and irregular event times
    data_continuous <- simulate_continuous_time_data(
        n = 1000,
        no_competing_events = FALSE,
        uncensored = FALSE
    )

    
    prep_data <- prepare_data(
        data = data_continuous,
        max_time_horizon = 720,
        time_covariates = c("A", "L"),
        baseline_covariates = c("age", "A_0", "L_0"),
        marginal_censoring = TRUE
    )
    altered_data <- propensity_scores(
        prepared_data = prep_data,
        model_treatment = "learn_glm_logistic",
        model_hazard = "learn_coxph"
    )

    # Run debiased ICE-IPCW procedure
    result <- debias_ice_ipcw(
        prepared_data = altered_data,
        time_horizon = 720,
        model_pseudo_outcome = "lm",
        model_hazard = "learn_coxph",
        conservative = TRUE,
        verbose = FALSE,
        penalize_pseudo_outcome = TRUE
    )

    correct_result <- data.table::data.table(
  estimate = 0.28487619124871394,
  se = 0.01684739758956828,
  lower = 0.25185529197316014,
  upper = 0.31789709052426773,
  ice_ipcw_estimate = 0.28585935059573064,
  ipw = 0.28495738884128674
)
    expect_true(all.equal(result, correct_result, tolerance = 1e-8))
})

test_that("test continuous time function (censored; competing events; ipcw_glm_expit, penalize)", {
    library(survival)
    library(data.table)
    library(prodlim)
    library(riskRegression)

    set.seed(34)
    # Simulate continuous time data with continuous and irregular event times
    data_continuous <- simulate_continuous_time_data(
        n = 1000,
        no_competing_events = FALSE,
        uncensored = FALSE
    )

    
    prep_data <- prepare_data(
        data = data_continuous,
        max_time_horizon = 720,
        time_covariates = c("A", "L"),
        baseline_covariates = c("age", "A_0", "L_0"),
        marginal_censoring = TRUE
    )
    altered_data <- propensity_scores(
        prepared_data = prep_data,
        model_treatment = "learn_glm_logistic",
        model_hazard = "learn_coxph"
    )

    # Run debiased ICE-IPCW procedure
    set.seed(65)
    result <- debias_ice_ipcw(
        prepared_data = altered_data,
        time_horizon = 720,
        model_pseudo_outcome = "ipcw_glm_expit",
        model_hazard = "learn_coxph",
        conservative = TRUE,
        verbose = FALSE,
        penalize_pseudo_outcome = TRUE
    )

    correct_result <- data.table::data.table(
  estimate = 0.2849064860164079,
  se = 0.016847871887799126,
  lower = 0.2518846571163216,
  upper = 0.3179283149164942,
  ice_ipcw_estimate = 0.28577586012994194,
  ipw = 0.28495738884128674
)
    expect_true(all.equal(result, correct_result, tolerance = 1e-8))
})


test_that("test continuous time function (censored; competing events; ipcw_glm_expit, penalize_treatment)", {
    library(survival)
    library(data.table)
    library(prodlim)
    library(riskRegression)

    set.seed(34)
    # Simulate continuous time data with continuous and irregular event times
    data_continuous <- simulate_continuous_time_data(
        n = 1000,
        no_competing_events = FALSE,
        uncensored = FALSE
    )

    
    prep_data <- prepare_data(
        data = data_continuous,
        max_time_horizon = 720,
        time_covariates = c("A", "L"),
        baseline_covariates = c("age", "A_0", "L_0"),
        marginal_censoring = TRUE
    )
    set.seed(65)
    altered_data <- propensity_scores(
        prepared_data = prep_data,
        model_treatment = "learn_glm_logistic",
        model_hazard = "learn_coxph",
        penalize_treatment = TRUE
    )

    # Run debiased ICE-IPCW procedure
    result <- debias_ice_ipcw(
        prepared_data = altered_data,
        time_horizon = 720,
        model_pseudo_outcome = "ipcw_glm_expit",
        model_hazard = "learn_coxph",
        conservative = TRUE,
        verbose = FALSE
    )

    correct_result <- data.table::data.table(
  estimate = 0.2847965969206015,
  se = 0.016833646114405392,
  lower = 0.2518026505363669,
  upper = 0.3177905433048361,
  ice_ipcw_estimate = 0.28426232748420016,
  ipw = 0.2852488617655428
)
    expect_true(all.equal(result, correct_result, tolerance = 1e-8))
})


test_that("test continuous time function (censored; competing events; ipcw_glm_expit, penalize_censoring)", {
    library(survival)
    library(data.table)
    library(prodlim)
    library(riskRegression)

    set.seed(34)
    # Simulate continuous time data with continuous and irregular event times
    data_continuous <- simulate_continuous_time_data(
        n = 1000,
        no_competing_events = FALSE,
        uncensored = FALSE
    )

    set.seed(65)
    prep_data <- prepare_data(
        data = data_continuous,
        max_time_horizon = 720,
        time_covariates = c("A", "L"),
        baseline_covariates = c("age", "A_0", "L_0"),
        marginal_censoring = FALSE
    )
    altered_data <- propensity_scores(
        prepared_data = prep_data,
        model_treatment = "learn_glm_logistic",
        model_hazard = "learn_coxph",
        penalize_hazard = TRUE
    )

    # Run debiased ICE-IPCW procedure
    result <- debias_ice_ipcw(
        prepared_data = altered_data,
        time_horizon = 720,
        model_pseudo_outcome = "oipcw_expit",
        model_hazard = "learn_coxph",
        conservative = TRUE,
        verbose = FALSE
    )

    correct_result <- data.table::data.table(
  estimate = 0.28474220931806654,
  se = 0.016843244444304774,
  lower = 0.2517294502072292,
  upper = 0.3177549684289039,
  ice_ipcw_estimate = 0.2846316598183721,
  ipw = 0.284940461612925
)
    expect_true(all.equal(result, correct_result, tolerance = 1e-8))
})
