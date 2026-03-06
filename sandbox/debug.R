devtools::load_all()
tar_source("~/phd/continuous_time_LTMLE/simulation_study/functions")
set.seed(14)

effects <- vary_effect(effect_A_on_Y =-0.3, effect_L_on_Y = 0.5, effect_L_on_A = -0.2, 
                         effect_A_on_L = -0.2, effect_age_on_Y = 0.025, effect_age_on_A = 0.002)
baseline_rate_list <- list(
    A = 0.005,
    L = 0.001,
    C = 0.0008,
    Y = 0.0001,
    D = 0.00015
)
K <- 3
data_continuous_intervention <- simulate_continuous_time_data(n = 100000, K = K, uncensored = FALSE, effects = effects, baseline_rate_list = baseline_rate_list, static_intervention = 1)$timevarying_data
data_continuous_intervention <- data_continuous_intervention[event != "A" & event != "L"]
true_val <- mean(1*(data_continuous_intervention$time <= 720 & data_continuous_intervention$event == "Y"))
res <- list()

apply_function <- function(args_debias_ice_ipcw, x, conservative = TRUE) {
        print(paste0("Runing with model_pseudo_outcome = ", x, " and conservative = ", conservative))
        if (x == "lm_penalized") {
            model_pseudo_outcome <- "lm"
            args_debias_ice_ipcw$penalize_pseudo_outcome <- TRUE
        } else {
            model_pseudo_outcome <- x
            args_debias_ice_ipcw$penalize_pseudo_outcome <- FALSE
        }
        args_debias_ice_ipcw$model_pseudo_outcome <- model_pseudo_outcome
        args_debias_ice_ipcw$conservative <- conservative
        res<-do.call(debias_ice_ipcw, args_debias_ice_ipcw)
        res$model_pseudo_outcome <- x
        res$conservative <- conservative
        res
    }

pseudo_outcomes <- c("lm","oipcw_expit", "ipcw_glm_expit", "lm_penalized") #"nls_expit", "nls_probit" "scaled_quasibinomial"
for (i in 1:1000) {
    temp <- list()
    set.seed(i)
    print(paste0("Iteration ", i))
    data_continuous_ate <- simulate_continuous_time_data(
        n = 1500,
        uncensored = FALSE,
        effects = effects,
        baseline_rate_list = baseline_rate_list,
        K = K
    )
    prep_data <- prepare_data(
        data = data_continuous_ate,
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
    args_debias_ice_ipcw <- list(
        prepared_data = altered_data,
        time_horizon = 720,
        model_hazard = "learn_coxph",
        conservative = FALSE,
        grid_size = 15
    )
    ## temp <- lapply(pseudo_outcomes, function(x) apply_function(args_debias_ice_ipcw, x, conservative = FALSE))
    ## temp <- rbindlist(temp)
    temp<-NULL
    args_debias_ice_ipcw$conservative <- TRUE
    temp_conservative <- lapply(pseudo_outcomes, function(x) apply_function(args_debias_ice_ipcw, x, conservative = TRUE))
    temp_conservative <- rbindlist(temp_conservative)
    temp2 <- data.table(
        estimate = true_val,
        se = NA,
        lower = NA,
        upper = NA,
        ice_ipcw_estimate = true_val,
        ipw = true_val,
        model_pseudo_outcome = "true_value",
        conservative = FALSE
    )
    res_temp <- rbindlist(list(temp2, temp, temp_conservative))
    setkeyv(res_temp, c("model_pseudo_outcome", "conservative"))
    res[[i]] <- res_temp
}
