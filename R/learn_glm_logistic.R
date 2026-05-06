### learn_glm_logistic.R --- 
#----------------------------------------------------------------------
## Author: Johan Sebastian Ohlendorff
## Created: Mar 13 2026 (18:55) 
## Version: 
## Last-Updated: May  6 2026 (11:52) 
##           By: Johan Sebastian Ohlendorff
##     Update #: 31
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

learn_glm_logistic <- function(character_formula, data, penalize, save_coefficients = FALSE, ...) {
  form <- as.formula(character_formula)
  coefficients <- NULL
  if (!isTRUE(penalize) || length(labels(stats::terms(form))) == 1L) {
    fit <- stats::glm(form, data = data, family = binomial(link = "logit"))
    if (save_coefficients) {
      coefficients <- coef(fit)
      coefficients[is.na(coefficients)] <- 0
    }
    predictions <- as.vector(predict(fit, type = "response"))
  } else {
    X <- model.matrix(form, data = data)[, -1]
    y <- data[[as.character(form[[2]])]]
    cv <- glmnet::cv.glmnet(X, y, alpha = 1, family = "binomial")
    fit <- glmnet::glmnet(X, y, alpha = 1, family = "binomial", lambda = cv$lambda.min)
    if (save_coefficients) {
      coefficients <- as.vector(coef(fit))
      names(coefficients) <- rownames(coef(fit)) 
    }
    predictions <- as.vector(predict(fit, newx = X, s = cv$lambda.min, type = "response"))
  }
  return(list(predictions = predictions, coefficients = coefficients))
}

######################################################################
### learn_glm_logistic.R ends here
