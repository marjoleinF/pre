context("Tests the pre functions including S3 methods")

test_that("Importance gives previous results with airquality data",{
  set.seed(42)
  airq.ens <- pre(Ozone ~ ., data=airquality[complete.cases(airquality),])
  
  imp1 <- importance(airq.ens, plot = FALSE)
  # save_to_test(imp1, "airquality_w_importance1")
  expect_equal(imp1, read_to_test("airquality_w_importance1"), tolerance = 1.490116e-08)
  
  imp2 <- importance(airq.ens, global = FALSE, plot = FALSE)
  # save_to_test(imp2, "airquality_w_importance2")
  expect_equal(imp2, read_to_test("airquality_w_importance2"), tolerance = 1.490116e-08)
  
  imp3 <- importance(airq.ens, global = FALSE, quantprobs = c(0, .25), plot = FALSE)
  # save_to_test(imp3, "airquality_w_importance3")
  expect_equal(imp3, read_to_test("airquality_w_importance3"), tolerance = 1.490116e-08)
})

test_that("Coef gives previous results with airquality data", {
  set.seed(42)
  airq.ens <- pre(Ozone ~ ., data=airquality[complete.cases(airquality),])
  
  coefs <- coef(airq.ens)
  
  # save_to_test(coefs, "airquality_w_pre_coef")
  expect_equal(coefs, read_to_test("airquality_w_pre_coef"), tolerance = 1.490116e-08)
})

test_that("cvpre gives previous results with airquality data", {
  set.seed(42)
  airq.ens <- pre(Ozone ~ ., data=airquality[complete.cases(airquality),])
  
  set.seed(7385056)
  airq.cv <- cvpre(airq.ens, k = 2)
  
  # save_to_test(airq.cv, "airquality_w_pre_cv")
  expect_equal(airq.cv, read_to_test("airquality_w_pre_cv"), tolerance = 1.490116e-08)
})

test_that("bsnullinteract and interact gives previous results with airquality data", {
  set.seed(42)
  airq.ens <- pre(Ozone ~ ., data=airquality[complete.cases(airquality),])
  
  set.seed(8969591)
  nullmods <- bsnullinteract(airq.ens, nsamp = 1)
  inter <- interact(airq.ens, c("Temp", "Wind"), nullmods = nullmods, plot = FALSE)
  
  for(i in 1:1)
    nullmods[[i]] <- nullmods[[i]][!names(nullmods[[i]]) %in%  c(
      "classify", "formula", "orig_data", "modmat_formula", "modmat", "data")]
  # save_to_test(nullmods, "airquality_w_bsnullinteract")
  expect_equal(nullmods, read_to_test("airquality_w_bsnullinteract"), tolerance = 1.490116e-08)
  
  # save_to_test(inter, "airquality_w_inter")
  expect_equal(inter, read_to_test("airquality_w_inter"), tolerance = 1.490116e-08)
})

test_that("Print gives previous resutls with airquality",{
  old <- getOption("digits")
  on.exit(options(digits = old))
  options(digits = 4)
  set.seed(42)
  airq.ens <- pre(Ozone ~ ., data=airquality[complete.cases(airquality),])
  
  # save_to_test(capture.output(print(airq.ens)), "airquality_print_output")
  expect_equal(capture.output(print(airq.ens)), read_to_test("airquality_print_output"), tolerance = 1.490116e-08)
  
  #####
  # Change digits
  options(digits = 2)
  
  # save_to_test(capture.output(print(airq.ens)), "airquality_print_output_less_digits")
  expect_equal(capture.output(print(airq.ens)), read_to_test("airquality_print_output_less_digits"), tolerance = 1.490116e-08)
})

test_that("Predict gives previous results with airquality data", {
  set.seed(42)
  airq.ens <- pre(Ozone ~ ., data=airquality[complete.cases(airquality),])
  
  preds <- predict(airq.ens, airquality[complete.cases(airquality),])
  # save_to_test(preds, "airquality_preds")
  
  expect_equal(preds, read_to_test("airquality_preds"), tolerance = 1.490116e-08)
})