test_that("gpre works with default settings and gives previous results", {
  #####
  # Continous outcome
  set.seed(8782650)
  fit <- gpre(
    Ozone ~ ., 
    data=airquality[complete.cases(airquality),])
  
  fit <- fit$glmnet.fit$glmnet.fit
  fit <- fit[c("a0", "beta")]
  fit$beta <- fit$beta[1:100, ]
  
  # save_to_test(fit, "gpre_fit1")
  expect_equal(fit, read_to_test("gpre_fit1"))
  
  #####
  # Binary outcome
  data(PimaIndiansDiabetes, package = "mlbench")
  
  set.seed(8782650)
  fit <- gpre(diabetes ~ ., PimaIndiansDiabetes, 
              base_learners = list(gpre_linear(), gpre_tress(ntrees = 20)))

  fit <- fit$glmnet.fit$glmnet.fit
  fit <- fit[c("a0", "beta")]
  fit$beta <- fit$beta[1:50, ]
  
  # save_to_test(fit, "gpre_fit1_binary")
  expect_equal(fit, read_to_test("gpre_fit1_binary"))
})

test_that("Sampling and subsampling works and is used in gpre", {
  #####
  # Bootstrap sample
  set.seed(seed <- 9495688)
  func <- gpre_sample(1)
  out <- func(100, rep(1, 100))
  
  expect_true(any(!1:100 %in% out))
  expect_length(out, 100)
  
  #####
  # Bootstrap with lower fraction
  func <- gpre_sample(.5)
  out <- func(100, rep(1, 100))
  
  expect_true(any(!1:100 %in% out))
  expect_true(any(table(out) > 1))
  expect_length(out, 50)
  
  #####
  # Bootstrap with weights
  func <- gpre_sample(1)
  ws <- rep(1, 100)
  ws[1] <- 100
  
  out <- func(100, ws)
  expect_true(names(sort(table(out),decreasing=TRUE))[1] == "1")
  expect_length(out, 100)
  
  #####
  # Bootstrap with lower fraction
  func <- gpre_sample(.5)
  out <- func(100, ws)
  
  expect_true(names(sort(table(out),decreasing=TRUE))[1] == "1")
  expect_length(out, 50)
})

test_that("gpre_tress gives expected_result for continous outcomes", {
  #####
  # Default settings with fewer trees
  func <- gpre_tress(ntrees = 100)
  
  args <- list(
    formula = Ozone ~ Solar.R + Wind + 
      cut(Wind, breaks = 3), # Notice that we include a factor. There where 
                             # some issues with factor at one point 
    data = airquality[complete.cases(airquality),],
    weights = rep(1, sum(complete.cases(airquality))),
    sample_func = gpre_sample(.5), 
    family = "gaussian")
  
  set.seed(seed <- 6772769)
  out <- do.call(func, args)
  
  # save_to_test(out, "gpre_tree_1")
  expect_equal(out, read_to_test("gpre_tree_1"), tolerance = 1.490116e-08)

  #####  
  # Only use stumps
  func <- gpre_tress(ntrees = 100, maxdepth = 1)
  
  set.seed(seed)
  out2 <- do.call(func, args)
  
  expect_lt(length(out2), length(out))
  expect_true(!any(grepl("\\$", out2)))
  # save_to_test(out2, "gpre_tree_2")
  expect_equal(out2, read_to_test("gpre_tree_2"), tolerance = 1.490116e-08)
  
  #####
  # Without learning rate
  func <- gpre_tress(ntrees = 100, learnrate = 0)
  
  set.seed(seed)
  out3 <- do.call(func, args)
  
  expect_true(any(!out3 %in% out))
  # save_to_test(out3, "gpre_tree_3")
  expect_equal(out3, read_to_test("gpre_tree_3"), tolerance = 1.490116e-08)
  
  #####
  # Without removal of duplicates and complements
  func <- gpre_tress(ntrees = 100, remove_duplicates_complements = FALSE)
  
  set.seed(seed)
  out4 <- do.call(func, args)
  
  expect_gt(length(out4), length(out))
  expect_true(all(out %in% out4))
  
  # save_to_test(out4, "gpre_tree_4")
  expect_equal(out4, read_to_test("gpre_tree_4"), tolerance = 1.490116e-08)
})

test_that("gpre_tress gives expected_result for binary outcomes", {
  #####
  # Works with binomial outcome
  data(PimaIndiansDiabetes, package = "mlbench")
  
  func <- gpre_tress(ntrees = 20)
  
  args <- list(
    formula = diabetes ~ ., 
    data = PimaIndiansDiabetes,
    weights = rep(1, nrow(PimaIndiansDiabetes)),
    sample_func = gpre_sample(), 
    family = "binomial")
  
  set.seed(seed <- 8779606)
  out <- do.call(func, args)
  
  # save_to_test(out, "gpre_tree_binary_1")
  expect_equal(out, read_to_test("gpre_tree_binary_1"), tolerance = 1.490116e-08)
  
  #####
  # Binary without learning rate
  func <- gpre_tress(ntrees = 100, learnrate = 0)
  
  set.seed(seed)
  out <- do.call(func, args)
  
  # save_to_test(out, "gpre_tree_binary_2")
  expect_equal(out, read_to_test("gpre_tree_binary_2"), tolerance = 1.490116e-08)
})

test_that("gpre_earth gives expected_result for continous outcomes", {
  ##### 
  # With defaults settings
  func <- gpre_earth()
  
  args <- list(
    formula = Ozone ~ Solar.R + Wind + 
      cut(Wind, breaks = 5), # Notice that we include a factor. There where  
                             # some issues with factor at one point 
    data = airquality[complete.cases(airquality),],
    weights = rep(1, sum(complete.cases(airquality))),
    sample_func = gpre_sample(.5), 
    family = "gaussian")
  
  set.seed(seed <- 6817505)
  out <- do.call(func, args)
  
  # save_to_test(out, "gpre_earth_1")
  expect_equal(out, read_to_test("gpre_earth_1"), tolerance = 1.490116e-08)
  
  ##### 
  # Additive model
  func <- gpre_earth(degree = 1)
  
  set.seed(seed)
  out2 <- do.call(func, args)
  
  # Match any string not containing *
  expect_match(out2, "^[^\\*]+$", perl = TRUE)
  # save_to_test(out2, "gpre_earth_2")
  expect_equal(out2, read_to_test("gpre_earth_2"), tolerance = 1.490116e-08)
  
  ##### 
  # Without high learning rate
  func <- gpre_earth(learnrate = 1)
  
  set.seed(seed)
  out3 <- do.call(func, args)
  
  expect_true(length(setdiff(out3, out)) > 0)
  # save_to_test(out3, "gpre_earth_3")
  expect_equal(out2, read_to_test("gpre_earth_2"), tolerance = 1.490116e-08)
  
  #####
  # With different number of end nodes
  func <- gpre_earth(nk = 5, ntrain = 1)
  
  set.seed(seed)
  out4 <- do.call(func, args)
  expect_length(out4, 4) # 5 - 1 for intercept
  
  func <- gpre_earth(nk = 50, ntrain = 1)
  set.seed(seed)
  out5 <- do.call(func, args)
  expect_gt(length(out5), length(out4))
  
  #####
  # Without standardizing
  func <- gpre_earth(standardize = FALSE)
  
  set.seed(seed)
  out6 <- do.call(func, args)
  
  expect_true(!any(grepl("scale\\ ?=", out6, perl = TRUE)))
  # save_to_test(out6, "gpre_earth_6")
  expect_equal(out6, read_to_test("gpre_earth_6"), tolerance = 1.490116e-08)
  
  ######
  # Continous outcome with defaults with two factors
  func <- gpre_earth(ntrain = 10)
  
  args <- list(
    formula = Ozone ~ Solar.R + Wind + 
      cut(Wind, breaks = 5) + # Notice that we include a factor. There where  
      cut(Temp, breaks = 5),  # some issues with factor at one point 
                                
    data = airquality[complete.cases(airquality),],
    weights = rep(1, sum(complete.cases(airquality))),
    sample_func = gpre_sample(.5), 
    family = "gaussian")
  
  set.seed(seed)
  out <- do.call(func, args)
})

test_that("gpre_earth gives expected_result for binary outcomes", {
  data(PimaIndiansDiabetes, package = "mlbench")
  
  func <- gpre_earth(ntrain = 10)
  
  args <- list(
    formula = diabetes ~ ., 
    data = PimaIndiansDiabetes,
    weights = rep(1, nrow(PimaIndiansDiabetes)),
    sample_func = gpre_sample(), 
    family = "binomial")
  
  set.seed(1067229)
  expect_message(out <- do.call(func, args), 
                 "Beware that gpre_earth will use L2 loss to train")
  
  # save_to_test(out, "gpre_earth_binary")
  expect_equal(out, read_to_test("gpre_earth_binary"), tolerance = 1.490116e-08)
})

test_that("gpre_linear gives expected_result", {
  func <- gpre_linear()
  
  args <- list(
    formula = Ozone ~ Solar.R + Wind + cut(Wind, breaks = 3), 
    data = airquality[complete.cases(airquality),],
    weights = rep(1, sum(complete.cases(airquality))))
  
  out <- do.call(func, args)
  
  # save_to_test(out, "gpre_linear_1")
  expect_equal(out, read_to_test("gpre_linear_1"), tolerance = 1.490116e-08)
  
  # No winsorization
  func <- gpre_linear(winsfrac = 0)
  
  out <- do.call(func, args)
  
  # save_to_test(out, "gpre_linear_2")
  expect_equal(out, read_to_test("gpre_linear_2"), tolerance = 1.490116e-08)
  
  # No scaling
  func <- gpre_linear(normalize = F)
  
  out <- do.call(func, args)
  
  # save_to_test(out, "gpre_linear_3")
  expect_equal(out, read_to_test("gpre_linear_3"), tolerance = 1.490116e-08)
  
  # do not do either
  func <- gpre_linear(normalize = F, winsfrac = 0)
  
  out <- do.call(func, args)
  
  expect_equal(out, c("lTerm(Solar.R)", "lTerm(Wind)"))
})