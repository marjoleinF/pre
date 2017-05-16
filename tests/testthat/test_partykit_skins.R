context("Tests the 'skin' versions of the partykit functions")

ctree_setup <- with(environment(pre), ctree_setup)
ctree_minmal <- with(environment(pre), ctree_minmal)

test_that("wrapper for ctree gives the same",{
  dat <- airquality[complete.cases(airquality),]
  
  set.seed(seed <- 2649390)
  org <- ctree(Ozone ~ ., data = dat, maxdepth = 3, mtry = Inf)
  
  set.seed(seed)
  tmp <- ctree_setup(Ozone ~ ., data = dat, maxdepth = 3, mtry = Inf)
  mini <- ctree_minmal(tmp$dat, tmp$response, tmp$control, tmp$ytrafo, tmp$terms)
  
  org <- unclass(org)
  mini <- unclass(mini)
  for(i in seq_along(org))
    expect_equal(org$node, mini$node)
  
  # save_to_test(mini$node, "tree_nodes_from_minimal")
  expect_equal(mini$node, read_to_test("tree_nodes_from_minimal"), tolerance = 0.00000001490116)
})

list.rules <- with(environment(pre), list.rules)

test_that("list.rules gives the same as previously",{
  dat <- airquality[complete.cases(airquality),]
  
  set.seed(2649390)
  tmp <- ctree_setup(Ozone ~ ., data = dat, maxdepth = 3, mtry = Inf)
  mini <- ctree_minmal(tmp$dat, tmp$response, tmp$control, tmp$ytrafo, tmp$terms)
  
  rules <- list.rules(mini)
  
  # save_to_test(rules, "tree_nodes_rules")
  expect_equal(rules, read_to_test("tree_nodes_rules"), tolerance = 0.00000001490116)
})

