### learn_glm_logistic.R --- 
#----------------------------------------------------------------------
## Author: Johan Sebastian Ohlendorff
## Created: Mar 13 2026 (18:55) 
## Version: 
## Last-Updated: Mar 25 2026 (14:29) 
##           By: Johan Sebastian Ohlendorff
##     Update #: 4
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

learn_glm_logistic <- function(character_formula, data, penalize, ...) {
  form <- as.formula(character_formula)
  if (!isTRUE(penalize) || length(labels(stats::terms(form))) == 1L) {
    fit <- stats::glm(form, data = data, family = binomial(link = "logit"))
    return(predict(fit, type = "response"))
  } else {
    X <- model.matrix(form, data = data)
    y <- data[[as.character(form[[2]])]]
    cv <- glmnet::cv.glmnet(X, y, alpha = 1, family = "binomial")
    fit <- glmnet::glmnet(X, y, alpha = 1, family = "binomial", lambda = cv$lambda.min)
    as.vector(predict(fit, newx = X, s = cv$lambda.min, type = "response"))
  }
}

######################################################################
### learn_glm_logistic.R ends here
