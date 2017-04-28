test_that("Print gives previous resutls with airquality",{
  set.seed(42)
  airq.ens <- pre(Ozone ~ ., data=airquality[complete.cases(airquality),])
  
  # save_to_test(capture.output(print(airq.ens)), "airquality_print_output")
  expect_equal(capture.output(print(airq.ens)), read_to_test("airquality_print_output"), tolerance = 1.490116e-08)
})