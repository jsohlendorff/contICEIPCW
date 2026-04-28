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
        time_horizons = 720,
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
                                      ipw = c(0.282975855447292),
                                      time_horizon = 720
                                  )

    expect_true(all.equal(result, correct_result, tolerance = 1e-8))
})

test_that("test continuous time function (censored; conservative)", {
    library(survival)
    library(data.table)
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
        time_horizons = 720,
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
                                      ipw = c(0.269324346522728),
                                      time_horizon = 720
                                  )

    expect_true(all.equal(result, correct_result, tolerance = 1e-8))
})


test_that("test continuous time function (censored; conservative, intervene_all_sequential_regression)", {
    library(survival)
    library(data.table)
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
        time_horizons = 720,
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
        model_pseudo_outcome = "scaled_quasibinomial",
        model_hazard = "learn_coxph",
        conservative = TRUE,
        verbose = FALSE,
        intervene_all_sequential_regression = TRUE
    )

    correct_result <- data.table::data.table(
                                      estimate = 0.27042666855710473,
                                      se = 0.01677307756052882,
                                      lower = 0.23755143653846825,
                                      upper = 0.30330190057574125,
                                      ice_ipcw_estimate = 0.2715374676034445,
                                      ipw = 0.26932434652272763,
                                      time_horizon = 720
                                  )
    expect_true(all.equal(result, correct_result, tolerance = 1e-8))
})



test_that("test continuous time function (censored; conservative, stratify)", {
    library(survival)
    library(data.table)
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
        time_horizons = 720,
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
        model_pseudo_outcome = "scaled_quasibinomial",
        model_hazard = "learn_coxph",
        conservative = TRUE,
        verbose = FALSE,
        stratify_sequential_regression = TRUE
    )

    correct_result <- data.table::data.table(
                                      estimate = 0.27018370034367295,
                                      se = 0.016761809897687625,
                                      lower = 0.23733055294420521,
                                      upper = 0.30303684774314066,
                                      ice_ipcw_estimate = 0.2712122534140999,
                                      ipw = 0.26932434652272763,
                                      time_horizon = 720
                                  )
    expect_true(all.equal(result, correct_result, tolerance = 1e-8))
})


test_that("test continuous time function (censored; conservative; marginal_censoring_hazard)", {
    library(survival)
    library(data.table)
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
        time_horizons = 720,
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
                                      ipw = 0.2693019050719549,
                                      time_horizon = 720
                                  )
    expect_true(all.equal(result, correct_result, tolerance = 1e-8))
})

test_that("test continuous time function (censored; non_conservative; multiple ice)", {
    library(survival)
    library(data.table)
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
        time_horizons = 720,
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
                                      ipw = 0.2693019050719549,
                                        time_horizon = 720
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
            time_horizons = 720,
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
            time_horizons = 720,
            time_covariates = c("A", "L"),
            baseline_covariates = c("age", "A_0", "L_0"),
            marginal_censoring = TRUE
        ),
        "There are ties in event times for some ids. Please ensure that each id has unique event times"
    )
})

test_that("update_TMLE option + version", {
    library(data.table)

    set.seed(34)
    # Simulate continuous time data with continuous and irregular event times
    data_continuous <- simulate_continuous_time_data(
        n = 1000,
        no_competing_events = TRUE,
        uncensored = FALSE
    )

    prep_data <- prepare_data(
        data = data_continuous,
        time_horizons = 720,
        time_covariates = c("A", "L"),
        baseline_covariates = c("age", "A_0", "L_0"),
        verbose = TRUE
    )
    altered_data <- propensity_scores(
        prepared_data = prep_data,
        model_treatment = "learn_glm_logistic",
        model_hazard = "learn_coxph",
        verbose = TRUE
    )

    # Run debiased ICE-IPCW procedure
    result <- debias_ice_ipcw(
                        prepared_data = altered_data,
                        model_pseudo_outcome = "oipcw_expit",
                        model_hazard = "learn_coxph",
                        conservative = TRUE,
                        verbose = TRUE,
                        tmle_update = TRUE
    )

    correct_result <- data.table::data.table(
                                      estimate = 0.2707014384434423,
                                      se = 0.01676985289712445,
                                      lower = 0.23783252676507838,
                                      upper = 0.3035703501218062,
                                      ice_ipcw_estimate = NA,
                                      ipw = 0.2693019050719549,
                                        time_horizon = 720
                                  )
    expect_true(all.equal(result, correct_result, tolerance = 1e-8))
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
        time_horizons = 720,
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
                                      ipw = 0.2593272643831233,
                                        time_horizon = 720
                                  )

    expect_true(all.equal(result, correct_result, tolerance = 1e-8))
})

test_that("test continuous time function (censored; conservative; competing risks)", {
    library(survival)
    library(data.table)
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
        time_horizons = 720,
        time_covariates = c("A", "L"),
        baseline_covariates = c("age", "A_0", "L_0"),
        marginal_censoring = FALSE
    )
    altered_data <- suppressWarnings(propensity_scores(
        prepared_data = prep_data,
        model_treatment = "learn_glm_logistic",
        model_hazard = "learn_coxph"
    ))

    # Run debiased ICE-IPCW procedure
    result <- debias_ice_ipcw(
        prepared_data = altered_data,
        
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
                                      ipw = 0.2852953695578807,
                                        time_horizon = 720
                                  )

    expect_true(all.equal(result, correct_result, tolerance = 1e-8))
})

test_that("test continuous time function (censored; competing events; conservative; marginal_censoring)", {
    library(survival)
    library(data.table)
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
        time_horizons = 720,
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
                                      ipw = 0.28495738884128674,
                                      time_horizon = 720
                                  )
    expect_true(all.equal(result, correct_result, tolerance = 1e-8))
})

test_that("test continuous time function (censored; competing events; conservative; marginal_censoring, multiple_ice", {
    library(survival)
    library(data.table)
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
        time_horizons = 720,
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
                                      ipw = 0.28495738884128674,
                                      time_horizon = 720
                                  )
    expect_true(all.equal(result, correct_result, tolerance = 1e-8))    
})

test_that("test continuous time function (censored; competing events; oicpw_expit)", {
    library(survival)
    library(data.table)
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
        time_horizons = 720,
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
                                      ipw = 0.28495738884128674,
                                      time_horizon = 720
                                  )
    expect_true(all.equal(result, correct_result, tolerance = 1e-8))
})

test_that("test continuous time function (censored; competing events; oicpw_probit)", {
    library(survival)
    library(data.table)
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
        time_horizons = 720,
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
                                      ipw = 0.28495738884128674,
                                      time_horizon = 720
                                  )
    expect_true(all.equal(result, correct_result, tolerance = 1e-8))
})

test_that("test continuous time function (censored; competing events; nls_expit", {
    library(survival)
    library(data.table)
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
        time_horizons = 720,
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
                                      ipw = 0.28495738884128674,
                                      time_horizon = 720
                                  )
    expect_true(all.equal(result, correct_result, tolerance = 1e-8))
})

test_that("test continuous time function (censored; competing events; nls_probit)", {
    library(survival)
    library(data.table)
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
        time_horizons = 720,
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
                                      ipw = 0.28495738884128674,
                                      time_horizon = 720
                                  )
    expect_true(all.equal(result, correct_result, tolerance = 1e-8))
})

test_that("test continuous time function (censored; competing events; ipcw_glm_expit)", {
    library(survival)
    library(data.table)
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
        time_horizons = 720,
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
                                      ipw = 0.28495738884128674,
                                      time_horizon = 720
                                  )
    expect_true(all.equal(result, correct_result, tolerance = 1e-8))
})

test_that("test continuous time function (censored; competing events; lm, penalize)", {
    library(survival)
    library(data.table)
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
        time_horizons = 720,
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
                                      ipw = 0.28495738884128674,
                                      time_horizon = 720
                                  )
    expect_true(all.equal(result, correct_result, tolerance = 1e-8))
})

test_that("test continuous time function (censored; competing events; ipcw_glm_expit, penalize)", {
    library(survival)
    library(data.table)
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
        time_horizons = 720,
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
                                      ipw = 0.28495738884128674,
                                      time_horizon = 720
                                  )
    expect_true(all.equal(result, correct_result, tolerance = 1e-8))
})


test_that("test continuous time function (censored; competing events; ipcw_glm_expit, penalize_treatment)", {
    library(survival)
    library(data.table)
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
        time_horizons = 720,
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
                                      ipw = 0.2852488617655428,
                                      time_horizon = 720
                                  )
    expect_true(all.equal(result, correct_result, tolerance = 1e-8))
})


test_that("test continuous time function (censored; competing events; ipcw_glm_expit, penalize_censoring)", {
    library(survival)
    library(data.table)
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
        time_horizons = 720,
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
                                      ipw = 0.284940461612925,
                                      time_horizon = 720
                                  )
    expect_true(all.equal(result, correct_result, tolerance = 1e-8))
})

test_that("test continuous time function (exclude variable)", {
    library(survival)
    library(data.table)
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
        time_horizons = 720,
        time_covariates = c("A", "L"),
        baseline_covariates = c("age", "A_0", "L_0"),
        marginal_censoring = TRUE
    )
    altered_data <- propensity_scores(
        prepared_data = prep_data,
        model_treatment = "learn_glm_logistic",
        model_hazard = "learn_coxph",
        penalize_hazard = FALSE,
        verbose = FALSE,
        exclude_latest_covariate = "time"
    )

    # Run debiased ICE-IPCW procedure
    result <- debias_ice_ipcw(
        prepared_data = altered_data,
        model_pseudo_outcome = "oipcw_expit",
        model_hazard = "learn_coxph",
        conservative = TRUE,
        verbose = FALSE
    )

    correct_result <- data.table::data.table(
                                      estimate = 0.28445004193932516,
                                      se = 0.016750137304035392,
                                      lower = 0.2516197728234158,
                                      upper = 0.3172803110552345,
                                      ice_ipcw_estimate = 0.284662875273469,
                                      ipw = 0.28368580294914786,
                                      time_horizon = 720
                                  )
    expect_true(all.equal(result, correct_result, tolerance = 1e-8))
})

test_that("test continuous time function (multiple time_horizons)", {
    library(survival)
    library(data.table)
    library(riskRegression)

    set.seed(34)
    # Simulate continuous time data with continuous and irregular event times
    data_continuous <- simulate_continuous_time_data(
        n = 1000,
        no_competing_events = FALSE,
        uncensored = FALSE
    )

    run_debias_ice_ipcw <- function(data_continuous, time_horizons) {
        prep_data <- prepare_data(
            data = data_continuous,
            time_horizons = time_horizons,
            time_covariates = c("A", "L"),
            baseline_covariates = c("age", "A_0", "L_0"),
            marginal_censoring = TRUE
        )
        altered_data <- propensity_scores(
            prepared_data = prep_data,
            model_treatment = "learn_glm_logistic",
            model_hazard = "learn_coxph",
            penalize_hazard = FALSE,
            verbose = FALSE,
            exclude_latest_covariate = "time"
        )

        # Run debiased ICE-IPCW procedure
        result <- debias_ice_ipcw(
            prepared_data = altered_data,
            model_pseudo_outcome = "oipcw_expit",
            model_hazard = "learn_coxph",
            conservative = TRUE,
            verbose = FALSE
        )
        return(result)
    }
    multiple_time_horizons_result <- run_debias_ice_ipcw(data_continuous, c(40, 720))
    result_40 <- run_debias_ice_ipcw(data_continuous, 40)
    result_720 <- run_debias_ice_ipcw(data_continuous, 720)
    results <- rbind(result_40, result_720)

    ## Check that results are the same as when running separately
    expect_true(all.equal(results, multiple_time_horizons_result, tolerance = 1e-8))
})


test_that("test continuous time function (multiple time_horizons; comparisons)", {
    library(survival)
    library(data.table)
    library(riskRegression)

    set.seed(34)
    # Simulate continuous time data with continuous and irregular event times
    data_continuous <- simulate_continuous_time_data(
        n = 1000,
        no_competing_events = FALSE,
        uncensored = FALSE
    )

    run_debias_ice_ipcw <- function(dat, time_horizons) {
        prep_data <- prepare_data(
            data = dat,
            time_horizons = time_horizons,
            time_covariates = c("A", "L"),
            baseline_covariates = c("age", "A_0", "L_0"),
            marginal_censoring = TRUE
        )
        altered_data <- propensity_scores(
            prepared_data = prep_data,
            model_treatment = "learn_glm_logistic",
            model_hazard = "learn_coxph",
            penalize_hazard = FALSE,
            verbose = FALSE,
            exclude_latest_covariate = "time"
        )

        # Run debiased ICE-IPCW procedure
        result <- debias_ice_ipcw(
            prepared_data = altered_data,
            model_pseudo_outcome = "oipcw_expit",
            model_hazard = "learn_coxph",
            conservative = TRUE,
            verbose = FALSE,
            return_ic = TRUE
        )
        return(result)
    }
    multiple_time_horizons_result <- run_debias_ice_ipcw(data_continuous, c(40, 720))
    multiple_time_horizons_result$treatment_name <- "A"
    ## Mess with data
    data_continuous_alt <- copy(data_continuous)
    set.seed(6545)
    data_continuous_alt$timevarying_data[event_number>1, A := rbinom(.N, 1, 0.8)]

    multiple_time_horizons_result_alt <- run_debias_ice_ipcw(data_continuous_alt, c(40, 720))
    multiple_time_horizons_result_alt$treatment_name <- "B"
    res<-compare_to_reference(
        reference_group = "A",
        multiple_time_horizons_result,
        multiple_time_horizons_result_alt
    )
    correct_result <- data.table::data.table(
                                      treatment = "B",
                                      reference_group = "A",
                                      time_horizon = c(40, 720),
                                      estimate = c(0, -0.000313726373610157),
                                      se = c(0, 0.008466625663265349),
                                      lower = c(0, -0.016908312673610242),
                                      upper = c(0, 0.016280859926389928)
                                  )
    expect_true(all.equal(res, correct_result, tolerance = 1e-8))
})

test_that("test continuous time function (reduce colinearity time)", {
    library(survival)
    library(data.table)
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
        time_horizons = 720,
        time_covariates = c("A", "L"),
        baseline_covariates = c("age", "A_0", "L_0"),
        marginal_censoring = TRUE
    )
    altered_data <- propensity_scores(
        prepared_data = prep_data,
        model_treatment = "learn_glm_logistic",
        model_hazard = "learn_coxph",
        penalize_hazard = FALSE,
        verbose = FALSE,
        reduce_colinearity_time = TRUE
    )

    # Run debiased ICE-IPCW procedure
    result <- debias_ice_ipcw(
        prepared_data = altered_data,
        model_pseudo_outcome = "oipcw_expit",
        model_hazard = "learn_coxph",
        conservative = TRUE,
        verbose = FALSE,
        reduce_colinearity_time = TRUE
    )

    correct_result <- data.table::data.table(
                                      estimate = 0.2847639451171848,
                                      se = 0.01684741059880183,
                                      lower = 0.2517430203435332,
                                      upper = 0.3177848698908364,
                                      ice_ipcw_estimate = 0.28466287527346895,
                                      ipw = 0.28495738884128674,
                                      time_horizon = 720
                                  )
    expect_true(all.equal(result, correct_result, tolerance = 1e-8))
})
