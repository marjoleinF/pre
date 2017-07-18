context("Tests the 'skin' versions of the partykit functions")

ctree_setup <- with(environment(pre), ctree_setup)
ctree_minimal <- with(environment(pre), ctree_minimal)
list.rules <- with(environment(pre), list.rules)

test_that("wrapper for ctree gives the same",{
  dat <- airquality[complete.cases(airquality),]
  
  set.seed(seed <- 2649390)
  org <- ctree(Ozone ~ ., data = dat, maxdepth = 3, mtry = Inf)
  
  set.seed(seed)
  tmp <- ctree_setup(Ozone ~ ., data = dat, maxdepth = 3, mtry = Inf)
  mini <- ctree_minimal(tmp$dat, tmp$response, tmp$control, tmp$ytrafo, tmp$terms)
  
  org <- unclass(org)
  mini <- unclass(mini)
  for(i in seq_along(org))
    expect_equal(org$node, mini$node)
  
  # save_to_test(mini$node, "tree_nodes_from_minimal")
  expect_equal(mini$node, read_to_test("tree_nodes_from_minimal"), tolerance = 0.00000001490116)
  
  #####
  # Had an error which motivated this test
  set.seed(6772769)
  tmp_dat <- airquality[complete.cases(airquality),]
  n <- nrow(tmp_dat)
  pick <- gpe_sample(.5)(n, rep(1, n))
  tmp_dat <- tmp_dat[pick, ]
  
  tmp <- ctree_setup(
    Ozone ~ Solar.R + Wind + cut(Wind, breaks = 3), 
    data = tmp_dat, maxdepth = 3, mtry = Inf)
  mini <- ctree_minimal(tmp$dat, tmp$response, tmp$control, tmp$ytrafo, tmp$terms)
  
  # save_to_test(unclass(mini), "tree_nodes_from_minimal_xtra")
  expect_equal(unclass(mini), read_to_test("tree_nodes_from_minimal_xtra"), tolerance = 1.490116e-08)
  
  # save_to_test(list.rules(mini), "tree_nodes_from_minimal_xtra_unlist")
  expect_equal(list.rules(mini), read_to_test("tree_nodes_from_minimal_xtra_unlist"), tolerance = 1.490116e-08)
})

test_that("list.rules gives the same as previously",{
  dat <- airquality[complete.cases(airquality),]
  
  set.seed(2649390)
  tmp <- ctree_setup(Ozone ~ ., data = dat, maxdepth = 3, mtry = Inf)
  mini <- ctree_minimal(tmp$dat, tmp$response, tmp$control, tmp$ytrafo, tmp$terms)
  
  rules <- list.rules(mini)
  
  # save_to_test(rules, "tree_nodes_rules")
  expect_equal(rules, read_to_test("tree_nodes_rules"), tolerance = 0.00000001490116)
})

