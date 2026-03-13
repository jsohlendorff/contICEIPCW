### learn_glm_logistic.R --- 
#----------------------------------------------------------------------
## Author: Johan Sebastian Ohlendorff
## Created: Mar 13 2026 (18:55) 
## Version: 
## Last-Updated: Mar 13 2026 (18:55) 
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
learn_glm_logistic <- function(character_formula,
                               data,
                               penalize, ...) {
  formula_object <- as.formula(character_formula)
  if (!penalize || length(labels(stats::terms(formula_object))) == 1){ ## do not run penalized regression with one covariate only
      ## Fit the logistic regression model
      fit <- stats::glm(formula_object, data = data, family = binomial(link = "logit"))
      ## Predict on original data
      predict(fit, type = "response")
  } else {
      ## Use Lasso with glmnet
      X <- model.matrix(formula_object, data = data)
      y <- data[[as.character(formula_object[[2]])]]
      cv_fit <- glmnet::cv.glmnet(X, y, alpha = 1, family = "binomial")
      fit <- glmnet::glmnet(X, y, alpha = 1, lambda = cv_fit$lambda.min, family = "binomial")
      as.vector(predict(fit, newx = X, s = "lambda.min", type = "response"))
  }  
}

######################################################################
### learn_glm_logistic.R ends here
