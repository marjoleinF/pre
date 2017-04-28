test_that("Predict gives previous results with airquality data", {
  set.seed(42)
  airq.ens <- pre(Ozone ~ ., data=airquality[complete.cases(airquality),])
  
  preds <- predict(airq.ens, airquality[complete.cases(airquality),])
  # save_to_test(preds, "airquality_preds")
  
  expect_equal(preds, read_to_test("airquality_preds"), tolerance = 1.490116e-08)
})