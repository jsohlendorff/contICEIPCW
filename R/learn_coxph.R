### learn_glm_logistic.R --- 
#----------------------------------------------------------------------
## Author: Johan Sebastian Ohlendorff
## Created: Mar 13 2026 (18:50) 
## Version: 
## Last-Updated: Mar 13 2026 (18:50) 
##           By: Johan Sebastian Ohlendorff
##     Update #: 1
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
# coxph learner for censoring
learn_coxph <- function(character_formula,
                        data,
                        time_variable = "time",
                        penalize, ...){
    exp_lp <- surv <- hazard <- NULL
    formula_cox <- as.formula(character_formula)

    if (!penalize || length(labels(stats::terms(as.formula(character_formula)))) == 1){ ## do not run penalized regression with one covariate only
        ## Fit the Cox model
        fit <- coxph(formula_cox, data = data, x = TRUE)
    } else {
        ## use glmnet
        fit<-riskRegression::GLMnet(formula_cox, data = data, family = "cox", alpha = 1)
    }
    list(pred = exp(-cumulative_hazard_cox(fit, data, data, time_variable, NULL)$Lambda_minus), fit = fit)
}

######################################################################
### learn_glm_logistic.R ends here
