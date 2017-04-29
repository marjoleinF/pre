test_that("Get previous results with airquality and pre function", {
  set.seed(42)
  airq.ens <- pre(Ozone ~ ., data=airquality[complete.cases(airquality),])
  
  # We remove some of the data to decrease the size
  airq.ens <- airq.ens[!names(airq.ens) %in%  c(
    "classify", "formula", "orig_data", "modmat_formula", "modmat", "data")]
  airq.ens$glmnet.fit <- airq.ens$glmnet.fit["glmnet.fit"]
  # save_to_test(airq.ens, "airquality_w_pre")
  expect_equal(airq.ens, read_to_test("airquality_w_pre"), tolerance = 1.490116e-08)
})

test_that("Get previous results with PimaIndiansDiabetes and pre function", {
  set.seed(7385056)
  data(PimaIndiansDiabetes, package = "mlbench")
  
  #####
  # Without learning rate
  fit <- pre(diabetes ~ ., data = PimaIndiansDiabetes)
  
  # We remove some of the data to decrease the size
  fit <- fit[!names(fit) %in%  c(
    "classify", "formula", "orig_data", "modmat_formula", "modmat", "data")]
  fit$glmnet.fit <- fit$glmnet.fit["glmnet.fit"]
  # save_to_test(fit, "PimaIndiansDiabetes_w_pre_no_LR")
  expect_equal(fit, read_to_test("PimaIndiansDiabetes_w_pre_no_LR"), tolerance = 1.490116e-08)
  
  #####
  # With learning rate
  fit <- pre(diabetes ~ ., data = PimaIndiansDiabetes, learnrate = .01,
             ntrees = 20, maxdepth = 3)
  
  # We remove some of the data to decrease the size
  fit <- fit[!names(fit) %in%  c(
    "classify", "formula", "orig_data", "modmat_formula", "modmat", "data", "rulevars")]
  fit$glmnet.fit <- fit$glmnet.fit["glmnet.fit"]
  # save_to_test(fit, "PimaIndiansDiabetes_w_pre_LR")
  expect_equal(fit, read_to_test("PimaIndiansDiabetes_w_pre_LR"), tolerance = 1.490116e-08)
})