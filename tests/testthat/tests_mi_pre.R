context("Tests multiple-imputation functions")


test_that("Get previous results with airquality and mi_pre function", {
  
  set.seed(42)
  
  ## Shoot holes in airquality data
  airq <- sapply(airquality, function(x) {
    x[sample(1:nrow(airquality), size = 25)] <- NA
    return(x)
  })
  
  ## Impute by replacing missing values with random sample from observed values
  imp <- replicate(3, data.frame(apply(airq, 2, function(x) {
    x[is.na(x)] <- sample(x[!is.na(x)], size = sum(is.na(x)))
    return(x)
    })), simplify = FALSE)
  
  ## Fit ensemble
  airq.ens <- mi_pre(Ozone ~ ., data = imp)
  # save_to_test(coef(airq.ens), "coef_airquality_w_mi_pre")
  expect_equal(coef(airq.ens), read_to_test("coef_airquality_w_mi_pre"), tolerance = 1.490116e-08)
})

test_that("Get previous results with airquality and mi_mean function", {
  set.seed(42)
  
  ## Shoot holes in airquality data
  airq <- sapply(airquality, function(x) {
    x[sample(1:nrow(airquality), size = 25)] <- NA
    return(x)
  })
  
  ## Impute by replacing missing values with random sample from observed values
  imp <- replicate(10, data.frame(apply(airq, 2, function(x) {
    x[is.na(x)] <- sample(x[!is.na(x)], size = sum(is.na(x)))
    return(x)
  })), simplify = FALSE)
  
  ## Create mean dataset
  mean_dat <- mi_mean(imp)
  # save_to_test(mean_dat, "airquality_w_mi_mean")
  expect_equal(mean_dat, read_to_test("airquality_w_mi_mean"), tolerance = 1.490116e-08)
})