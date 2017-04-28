test_that("Get previous results with airquality and pre function", {
  set.seed(42)
  airq.ens <- pre(Ozone ~ ., data=airquality[complete.cases(airquality),])
  
  # We remove some of the data to decrease the size
  airq.ens <- airq.ens[!names(airq.ens) %in%  c(
    "classify", "formula", "orig_data", "modmat_formula", "modmat", "data")]
  # save_to_test(airq.ens, "airquality_w_pre")
  expect_equal(airq.ens, read_to_test("airquality_w_pre"), tolerance = 1.490116e-08)
})