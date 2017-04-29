ctree_setup <- with(environment(pre), ctree_setup)
ctree_minmal <- with(environment(pre), ctree_minmal)

test_that("wrapper for ctree gives the same",{
  dat <- airquality[complete.cases(airquality),]
  
  set.seed(seed <- 2649390)
  org <- ctree(Ozone ~ ., data = dat, maxdepth = 3, mtry = Inf)
  
  set.seed(seed)
  tmp <- ctree_setup(Ozone ~ ., data = dat, maxdepth = 3, mtry = Inf)
  mini <- ctree_minmal(tmp$dat, tmp$response, tmp$control, tmp$ytrafo)
  
  org <- unclass(org)
  mini <- unclass(mini)
  for(i in seq_along(org))
    expect_equal(org$node, mini$node)
})
