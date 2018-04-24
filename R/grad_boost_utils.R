#####
# Functions for gradient boosting

get_intercept_logistic <- function(y, ws = NULL) {
  # # page 484 of:
  # # Bühlmann, Peter, and Torsten Hothorn. "Boosting algorithms: Regularization, 
  # # prediction and model fitting." Statistical Science (2007): 477-505.
  # # or check do the math an figure out that:
  # n <- 1000
  # y <- runif(n) > 1/(1 + exp(-1))
  # w <- runif(n, 0, 2)
  # 
  # glm.fit(
  #   matrix(rep(1, n), ncol = 1), 
  #   y,
  #   family = binomial(), 
  #   weights = w)$coefficients
  # 
  # p <- weighted.mean(y, w)
  # log(p / (1 - p))
  
  p_bar <- if(is.null(ws)) mean(y) else weighted.mean(y, ws)
  log(p_bar / (1 - p_bar))
}

get_y_learn_logistic <- function(eta, y) {
  # See LogitBoost on page 351 of:
  # Friedman, J., Hastie, T., & Tibshirani, R. (2000). Additive logistic 
  # regression: a statistical view of boosting (with discussion and a rejoinder 
  # by the authors). The annals of statistics, 28(2), 337-407.
  
  trunc_fac <- 12
  eta <- pmin(pmax(eta, -trunc_fac), trunc_fac)
  p <- 1 / (1 + exp(-eta))
  (y - p) / sqrt(p * (1 - p))
}

get_intercept_count <- function(y, ws = NULL) {
  lambda_bar <- if(is.null(ws)) mean(y) else weighted.mean(y, ws)
  log(lambda_bar)
}

get_y_learn_count <- function(eta, y) {
  lambda <- exp(eta)
  y - lambda
}


get_intercept_multinomial <- function(y, ws = NULL) {
  p_bar <- if (is.null(ws)) colMeans(y) else apply(y, 2, mean, weights = ws)
  log(p_bar / (1 - p_bar))
}

get_y_learn_multinomial <- function(eta, y) {
  #eta <- cbind(-14:15, 14:-15)
  #y <- cbind(rep(0:1, each = 15), rep(1:0, each = 15))
  eta <- apply(eta, 2, function(eta, trunc_fac = 12) {
    pmin(pmax(eta, -trunc_fac), trunc_fac)
  })
  p <- 1 / (1 + exp(-eta))
  (y - p) / sqrt(p * (1 - p))
}


