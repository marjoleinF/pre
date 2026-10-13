context("Tests functions explain and shap.")

test_that("functions explain and shap give previous results with factor and numeric predictors", {
  set.seed(42)
  airq <- airquality
  airq$Month <- factor(airq$Month)
  airq$Day <- factor(airq$Day)
  airq.ens <- pre(Ozone ~ ., data=airq, ntrees = 10)
  expl <- explain(airq.ens, newdata = airq, deprecation.message = FALSE)
  # save_to_test(expl, "explain")
  expect_equal(expl, read_to_test("explain"))

  SHAP <- shap(airq.ens, newdata = airq[1:10, ], reference_data = airq)  
  # save_to_test(SHAP$marginal, "SHAP_marginal")
  expect_equal(SHAP$marginal, read_to_test("SHAP_marginal"))
  # save_to_test(SHAP$interactions, "SHAP_interact")
  expect_equal(SHAP$interactions, read_to_test("SHAP_interact"))
})
