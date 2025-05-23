context("Tests the pre functions")


test_that("Get previous results with airquality and pre function", {
  set.seed(42)
  #####
  # With learning rate
  airq.ens <- pre(Ozone ~ ., data=airquality, ntrees = 10)
  
  # We remove some of the data to decrease the size
  airq.ens <- airq.ens[!names(airq.ens) %in%  c(
    "formula", "modmat_formula", "modmat", "data")]
  airq.ens$glmnet.fit <- airq.ens$glmnet.fit["glmnet.fit"]
  # save_to_test(airq.ens, "airquality_w_pre")
  expect_equal(airq.ens, read_to_test("airquality_w_pre"), tolerance = 1.490116e-08)
  
  #####
  # Without learning rate
  set.seed(42)
  airq.ens <- pre(Ozone ~ ., data=airquality, learnrate = 0, ntrees = 10)
  airq.ens <- airq.ens[!names(airq.ens) %in%  c(
    "formula", "modmat_formula", "modmat", "data")]
  airq.ens$glmnet.fit <- airq.ens$glmnet.fit["glmnet.fit"]
  # save_to_test(airq.ens, "airquality_w_pre_without_LR")
  expect_equal(airq.ens, read_to_test("airquality_w_pre_without_LR"), tolerance = 1.490116e-08)
  
  #####
  # Works only with rules
  set.seed(42)
  airq.ens <- pre(Ozone ~ ., data = airquality, type="rules", ntrees = 10)
  airq.ens <- airq.ens[!names(airq.ens) %in%  c(
    "formula", "modmat_formula", "modmat", "data")]
  airq.ens$glmnet.fit <- airq.ens$glmnet.fit["glmnet.fit"]
  # save_to_test(airq.ens, "airquality_w_pre_without_linear_terms")
  expect_equal(airq.ens, read_to_test("airquality_w_pre_without_linear_terms"), tolerance = 1.49e-08)
  
  ####
  # Works with multivariate responses
  set.seed(42)
  airq.ens <- pre(Ozone + Solar.R ~ ., data = airquality, family = "mgaussian", ntrees = 10)
  airq.ens <- airq.ens[!names(airq.ens) %in%  c(
    "formula", "modmat_formula", "modmat", "data")]
  airq.ens$glmnet.fit <- airq.ens$glmnet.fit["glmnet.fit"]
  # save_to_test(airq.ens, "airquality_w_pre_with_multivariate_response")
  expect_equal(airq.ens, read_to_test("airquality_w_pre_with_multivariate_response"), tolerance = 1.49e-08)

  
  #####
  # Works with rpart
  set.seed(42)
  airq.ens <- pre(Ozone ~ ., data = airquality, tree.unbiased = FALSE, ntrees = 10)
  # save_to_test(airq.ens$rules, "airquality_w_pre_with_rpart")
  expect_equal(airq.ens$rules, read_to_test("airquality_w_pre_with_rpart"), tolerance = 1.49e-08)
  
  #####
  # Works with glmtree
  set.seed(42)
  airq.ens <- pre(Ozone ~ ., data = airquality, use.grad = FALSE, ntrees = 10)
  # save_to_test(airq.ens, "airquality_w_pre_with_glmtree")
  expect_equal(airq.ens, read_to_test("airquality_w_pre_with_glmtree"), tolerance = 1.49e-08)
  
  ####
  # (hidden) argument offset is passed to function glmtree and cv.glmnet
  set.seed(42)
  airq.ens.offs <- pre(Ozone ~ ., data = airquality, use.grad = FALSE, ntrees = 10, 
                       offset = airquality$Wind)
  expect(nrow(airq.ens$rules) != nrow(airq.ens.offs$rules), 
         failure = "Use of offset does not affect the number of generated rules (which it previously did.")
  expect(is.null(airq.ens.offs$glmnet.fit$call$offset) != is.null(airq.ens$glmnet.fit$call$offset),
         failure = "Argument offset might not be correcltly passed to cv.glmnet")
  # save_to_test(airq.ens.offs, "airquality_w_pre_with_glmtree_w_offset")
  expect_equal(airq.ens.offs, read_to_test("airquality_w_pre_with_glmtree_w_offset"), tolerance = 1.49e-08)
  
  #####
  # Works with adaptive lasso
  set.seed(42)
  airq.ens <- pre(Ozone ~ ., data = airquality, ad.alpha = 0, ntrees = 10)
  airq.ens$glmnet.fit$call <- airq.ens$glmnet.fit$glmnet.fit$call <- NULL
  # save_to_test(airq.ens, "airquality_w_pre_with_adaptive_lasso")
  expect_equal(airq.ens, read_to_test("airquality_w_pre_with_adaptive_lasso"), tolerance = 1.49e-08)
  
  ####
  # Works with relaxed adaptive lasso
  set.seed(42)
  airq.ens <- pre(Ozone ~ ., data = airquality, ad.alpha = 0, relax = TRUE, ntrees = 10)
  airq.ens$glmnet.fit$call <- airq.ens$glmnet.fit$glmnet.fit$call <- NULL
  airq.ens$glmnet.fit$glmnet.fit$relaxed$call <- NULL
  # save_to_test(airq.ens, "airquality_w_pre_with_relaxed_adaptive_lasso")
  expect_equal(airq.ens, read_to_test("airquality_w_pre_with_relaxed_adaptive_lasso"), tolerance = 1.49e-08)
  
  #### 
  # Relaxed (adaptive) lasso works with function prune_pre
  a <- capture.output(prune_pre(airq.ens, nonzero = 2))
  # save_to_test(a, "airquality_w_relaxed_pre_with_prune_pre")
  expect_equal(a, read_to_test("airquality_w_relaxed_pre_with_prune_pre"), tolerance = 1.49e-08)
  
})


test_that("Get previous results with PimaIndiansDiabetes and pre function", {
  #####
  # Without learning rate
  set.seed(7385056)
  fit <- pre(diabetes ~ ., data = PimaIndiansDiabetes, ntrees = 10, maxdepth = 3)
  
  # We remove some of the data to decrease the size
  fit <- fit[names(fit) %in%  c("rules", "glmnet.fit")]
  fit$glmnet.fit <- fit$glmnet.fit["glmnet.fit"]
  fit$glmnet.fit$glmnet.fit <- fit$glmnet.fit$glmnet.fit["beta"]
  fit$glmnet.fit$glmnet.fit[["beta"]] <- head(fit$glmnet.fit$glmnet.fit$beta@x, 200)
  fit$rules <- as.matrix(fit$rules)
  # save_to_test(fit, "PimaIndiansDiabetes_w_pre_no_LR")
  expect_equal(fit, read_to_test("PimaIndiansDiabetes_w_pre_no_LR"), tolerance = 1.490116e-08)
  
  #####
  # With learning rate
  set.seed(4989935)
  fit <- pre(diabetes ~ ., data = PimaIndiansDiabetes, learnrate = .01,
             ntrees = 10, maxdepth = 3)
  
  # We remove some of the data to decrease the size
  fit <- fit[!names(fit) %in%  c(
    "formula", "modmat_formula", "modmat", "data")]
  fit$glmnet.fit <- fit$glmnet.fit["glmnet.fit"]
  # save_to_test(fit, "PimaIndiansDiabetes_w_pre_LR")
  expect_equal(fit, read_to_test("PimaIndiansDiabetes_w_pre_LR"), tolerance = 1.490116e-08)
})



test_that("Get previous results with iris and pre function", {
  
  #####
  # With learning rate
  set.seed(7385056)
  fit <- pre(Species ~ ., data = iris, ntrees = 10, maxdepth = 3)
  
  # We remove some of the data to decrease the size
  fit <- fit[names(fit) %in%  c("rules", "glmnet.fit")]
  fit$glmnet.fit <- fit$glmnet.fit["glmnet.fit"]
  fit$glmnet.fit$glmnet.fit <- fit$glmnet.fit$glmnet.fit["beta"]
  fit$glmnet.fit$glmnet.fit[["beta"]]$virginica <- fit$glmnet.fit$glmnet.fit[["beta"]]$virginica
  fit$rules <- as.matrix(fit$rules)
  # save_to_test(fit, "iris_w_pre")
  expect_equal(fit, read_to_test("iris_w_pre"), tolerance = 1.490116e-06)
  
  if (Sys.info()["sysname"] != "SunOS") {
    ##
    ## These tests fail on solaris only. Cannot test whether this still fails, so commented out for now.
    ## See https://cran.r-project.org/web/checks/check_results_pre.html
    ## Could help to use 0L instead of 0 for learning rate?
    ##
    #####
    # Without learning rate
    set.seed(4989935)
    fit <- pre(Species ~ ., data = iris, learnrate = 0, ntrees = 10, maxdepth = 3, nlambda =10)
    
    # We remove some of the data to decrease the size
    fit <- fit[names(fit) %in%  c("rules", "glmnet.fit")]
    fit$call <- NULL
    fit$glmnet.fit <- fit$glmnet.fit["glmnet.fit"]
    fit$glmnet.fit$glmnet.fit <- fit$glmnet.fit$glmnet.fit["beta"]
    fit$glmnet.fit$glmnet.fit[["beta"]]$virginica <- fit$glmnet.fit$glmnet.fit[["beta"]]$virginica
    fit$rules <- as.matrix(fit$rules)
    # save_to_test(fit, "iris_w_pre_no_learn")
    expect_equal(fit, read_to_test("iris_w_pre_no_learn"), tolerance = 1.490116e-06)
    
    #####
    # Without learning rate, parallel computation
    library("doParallel")
    cl <- makeCluster(2)
    registerDoParallel(cl)
    set.seed(4989935)
    fit2 <- pre(Species ~ ., data = iris, learnrate = 0, ntrees = 10, maxdepth = 3, par.init=TRUE, par.final=TRUE, nlambda = 10)
    stopCluster(cl)
    # We remove some of the data to decrease the size
    fit2 <- fit2[names(fit2) %in%  c("rules", "glmnet.fit")]
    fit2$call <- NULL
    fit2$glmnet.fit <- fit2$glmnet.fit["glmnet.fit"]
    fit2$glmnet.fit$glmnet.fit <- fit2$glmnet.fit$glmnet.fit["beta"]
    fit2$glmnet.fit$glmnet.fit[["beta"]]$virginica <- fit2$glmnet.fit$glmnet.fit[["beta"]]$virginica
    fit2$rules <- as.matrix(fit2$rules)
    # save_to_test(fit, "iris_w_pre_no_learn_par")
    expect_equal(fit, fit2)
    expect_equal(fit2, read_to_test("iris_w_pre_no_learn_par"), tolerance = 1.490116e-06)
  }
})

test_that("Get previous results with lung survival data", {
    set.seed(42)
    fit <- pre(Surv(time, status) ~ ., data = Lung, ntrees = 10, family = "cox")
    fit <- fit[names(fit) %in%  c("rules", "glmnet.fit")]
    fit$call <- NULL
    fit$glmnet.fit <- fit$glmnet.fit["glmnet.fit"]
    fit$glmnet.fit$glmnet.fit <- fit$glmnet.fit$glmnet.fit["beta"]
    fit$rules <- as.matrix(fit$rules)
    # save_to_test(fit, "lung_w_pre_surv")
    expect_equal(fit, read_to_test("lung_w_pre_surv"), tolerance = 1.490116e-06)
})

test_that("pre gives almost the same rules with `sparse` set to `TRUE` and `FALSE`", {
  #####
  # With learning rate
  set.seed(seed <- 42)
  airq.ens <- pre(Ozone ~ ., data=airquality, ntrees = 10)
  set.seed(seed)
  airq.ens.sparse <- pre(Ozone ~ ., data=airquality, ntrees = 10, sparse = TRUE)
  
  # will differ because sparse complements are used
  expect_false(isTRUE(all.equal(
    airq.ens$glmnet.fit$glmnet.fit$beta, 
    airq.ens.sparse$glmnet.fit$glmnet.fit$beta)))
  
  expect_false(all(
    airq.ens$rules$description == 
      airq.ens.sparse$rules$description))
  
  # but the prediction are identical 
  expect_equal(predict(airq.ens), predict(airq.ens.sparse))
  
  # check the classes of their model matrices
  expect_equal(inherits(airq.ens$modmat, "matrix"), TRUE)
  expect_s4_class(airq.ens.sparse$modmat, "dgCMatrix")
})

test_that("predict.pre works with `sparse` set to `TRUE`", {
  #####
  # With learning rate
  set.seed(42)
  airq.ens.sparse <- pre(Ozone ~ ., data=airquality, ntrees = 10, sparse = TRUE)
  
  expect_equal(
    predict(airq.ens.sparse, type = "link", newdata = airquality[1:10, ]),
    predict(airq.ens.sparse, type = "link")[1:10])
})

test_that("Get previous results with iris and rare_level_sampler function", {
  ## Create dataset with two factors containing rare levels
  dat <- iris[iris$Species != "versicolor", ]
  dat <- rbind(dat, iris[iris$Species == "versicolor", ][1:5, ])
  dat$factor2 <- factor(rep(1:21, times = 5))
  ## Obtain samples
  samp_func <- rare_level_sampler(c("Species", "factor2"), data = dat, 
                                  sampfrac = .2)
  N <- nrow(dat)
  wts <- rep(1, times = nrow(dat))
  set.seed(3)
  samples <- list()
  for (i in 1:10) samples[[i]] <- dat[samp_func(n = N, weights = wts), ]
  # save_to_test(samples, "rare_level_sampler")
  expect_equal(samples, read_to_test("rare_level_sampler"), tolerance = 1.490116e-08)  
})