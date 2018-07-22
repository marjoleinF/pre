context("Tests the pre functions with S3 methods and maxdepth sampler.")

test_that("maxdepth.sampler gives previous results", {
  
  func1 <- maxdepth_sampler()
  set.seed(42)
  maxdepths <- func1(100)
  # save_to_test(maxdepths, "maxdepths_default")
  expect_equal(maxdepths, read_to_test("maxdepths_default"), tolerance = 1.490116e-03)
  
  func2 <- maxdepth_sampler(av.no.term.nodes = 16L)
  set.seed(42)
  maxdepths2 <- func2(100)
  # save_to_test(maxdepths2, "maxdepths_16_term_nodes")
  expect_equal(maxdepths2, read_to_test("maxdepths_16_term_nodes"), tolerance = 1.490116e-03)
  
  func3 <- maxdepth_sampler(av.tree.depth = 7)
  set.seed(42)
  maxdepths3 <- func3(100)
  # save_to_test(maxdepths3, "maxdepths_tree_depth_7")
  expect_equal(maxdepths3, read_to_test("maxdepths_tree_depth_7"), tolerance = 1.490116e-03)
})


test_that("Importance gives previous results with airquality data",{
  set.seed(42)
  airq.ens <- pre(Ozone ~ ., data=airquality, ntrees = 10)
  
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
  airq.ens <- pre(Ozone ~ ., data=airquality, ntrees = 10)
  
  coefs <- coef(airq.ens)
  
  # save_to_test(coefs, "airquality_w_pre_coef")
  expect_equal(coefs, read_to_test("airquality_w_pre_coef"), tolerance = 1.490116e-08)
})

test_that("cvpre gives previous results with airquality data", {
  set.seed(42)
  airq.ens <- pre(Ozone ~ ., data=airquality, ntrees = 10)
  
  set.seed(7385056)
  airq.cv <- cvpre(airq.ens, k = 2)
  
  # save_to_test(airq.cv, "airquality_w_pre_cv")
  expect_equal(airq.cv, read_to_test("airquality_w_pre_cv"), tolerance = 1.490116e-08)
})

test_that("bsnullinteract and interact gives previous results with airquality data", {
  set.seed(42)
  airq.ens <- pre(Ozone ~ ., data=airquality, ntrees = 10)
  
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

test_that("Print gives previous results with airquality",{
  old <- getOption("digits")
  on.exit(options(digits = old))
  options(digits = 4)
  set.seed(42)
  airq.ens <- pre(Ozone ~ ., data=airquality, ntrees = 10)
  
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
  airq.ens <- pre(Ozone ~ ., data=airquality, ntrees = 10)
  preds <- predict(airq.ens, airquality)

  # plot(preds, airquality$Ozone)
  # abline(a = 0, b = 1)
  # mean((preds - airquality$Ozone)^2)
  
  # save_to_test(preds, "airquality_preds")
  expect_equal(preds, read_to_test("airquality_preds"), tolerance = 1.490116e-08)
  
  #####
  # Should give the same data to use in fit is not passed
  other_preds <- predict(airq.ens)
  
  expect_equal(other_preds, preds)
})

test_that("`delete_duplicates_complements` gives the same regardless of `sparse` argument", {
  set.seed(25106787)
  X <- data.frame(
    matrix(runif(10 * 20), 10, dimnames = list(NULL, paste0("X", 1:20))))
  rules <- paste0(
    "X", sample.int(20, 100, replace = TRUE), 
    ifelse(runif(100) > .5, ">", "<="), runif(100))
  
  confs <- expand.grid(
    removecomplements = c(T, F), removeduplicates = c(T, F), 
    return.dupl.compl = c(T, F))
  
  for(i in 1:nrow(confs)) { 
    v <- confs[i, ]
    eval(bquote(
      expect_equal(
        delete_duplicates_complements(
          removecomplements = .(v$removecomplements), 
          removeduplicates  = .(v$removeduplicates), 
          return.dupl.compl = .(v$return.dupl.compl),
          rules = rules, data = X, sparse = TRUE),
        delete_duplicates_complements(
          removecomplements = .(v$removecomplements), 
          removeduplicates  = .(v$removeduplicates), 
          return.dupl.compl = .(v$return.dupl.compl),
          rules = rules, data = X))))
  }
})