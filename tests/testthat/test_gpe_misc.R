context("Tests the gpe S3 functions")

test_that("Print works for gpe", {
  dig <- getOption("digits")
  on.exit(options(digits = dig))
  options(digits = 3)
  
  #####
  # W/ linear terms
  suppressWarnings(RNGversion("3.1.0"))
  set.seed(4825707)
  airq.ens <- gpe(
    Ozone ~ .,
    data= airquality,
    base_learners = list(gpe_linear()))
  
  # save_to_test(capture.output(print(airq.ens)), "gpe_airquality_print_output_linear")
  expect_equal(capture.output(print(airq.ens)), read_to_test("gpe_airquality_print_output_linear"))
  
  #####
  # W/ trees
  airq.ens <- gpe(
    Ozone ~ .,
    data= airquality,
    base_learners = list(gpe_trees(ntrees = 3)))
  
  # save_to_test(capture.output(print(airq.ens)), "gpe_airquality_print_output")
  expect_equal(capture.output(print(airq.ens)), 
               read_to_test("gpe_airquality_print_output"))
  
  #####
  # W/ earth
  airq.ens <- gpe(
    Ozone ~ .,
    data= airquality,
    base_learners = list(
      gpe_earth(
        ntrain = 3)))
  
  # save_to_test(capture.output(print(airq.ens)), "gpe_airquality_print_output_earth")
  expect_equal(capture.output(print(airq.ens)), read_to_test("gpe_airquality_print_output_earth"))
})

test_that("Coef works for gpe", {
  suppressWarnings(RNGversion("3.1.0"))
  set.seed(9116073)
  airq.ens <- gpe(
    Ozone ~ .,
    data=airquality,
    base_learners = list(gpe_trees(ntrees = 10)))
  
  coefs <- coef(airq.ens)
  
  # save_to_test(coefs, "gpe_airquality_w_pre_coef")
  expect_equal(coefs, read_to_test("gpe_airquality_w_pre_coef"), tolerance = 1.49e-08)
})

test_that("Predict works for gpe and gives previous results", {
  #####
  # Regression problem
  suppressWarnings(RNGversion("3.1.0"))
  set.seed(seed <- 9638602)
  airq.ens <- gpe(
    Ozone ~ .,
    data=airquality,
    base_learners = list(gpe_trees(ntrees = 10)))
  
  preds <- predict(airq.ens)
  # plot(preds ~ airquality$Ozone[complete.cases(airquality)])
  # abline(a = 0, b = 1, lty = 2)
  
  # save_to_test(preds, "gpe_predict_regression")
  expect_equal(preds, read_to_test("gpe_predict_regression"), tolerance = 1.49e-08)
  
  suppressWarnings(RNGversion("3.1.0"))
  set.seed(seed)
  airq.ens <- gpe(
    Ozone ~ .,
    data=airquality, 
    base_learners = list(gpe_trees(ntrees = 10)),
    model = FALSE) # don't save model
  
  expect_null(airq.ens$data)
  expect_error(predict(airq.ens), 
               "Predict called with no new object and no saved data with gpe")
  
  preds_new <- predict(
    airq.ens, newdata = airquality)
  expect_equal(preds, preds_new)
  
  #####
  # Binary
  suppressWarnings(RNGversion("3.1.0"))
  set.seed(seed)
  fit <- gpe(
    diabetes ~ ., data = PimaIndiansDiabetes,
    base_learners = list(gpe_trees(ntrees = 3), gpe_linear()))
  
  preds <- predict(fit, type = "response")
  
  # actual <- PimaIndiansDiabetes$diabetes == levels(PimaIndiansDiabetes$diabetes)[2]
  # plot(smooth.spline(preds, actual), type = "l", ylim = c(-.2, 1.2), col = "red")
  # points(preds, jitter(as.integer(actual), .1))
  # hist(preds)
  
  # save_to_test(preds[1:100], "gpe_predict_binary_response")
  expect_equal(preds[1:100], read_to_test("gpe_predict_binary_response"), tolerance = 1.49e-08)
})
