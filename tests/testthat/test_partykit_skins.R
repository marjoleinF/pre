context("Tests the 'skin' versions of the partykit functions")

test_that("list.rules gives the same as previously",{
  dat <- airquality[complete.cases(airquality),]
  
  set.seed(2649390)
  tree <- ctree(Ozone ~ ., data = dat, 
                control = ctree_control(maxdepth = 3, mtry = Inf))
  
  rules <- list.rules(tree)
  
  expect_length(rules, 2^3)
  
  # save_to_test(rules, "tree_nodes_rules")
  expect_equal(rules, read_to_test("tree_nodes_rules"), tolerance = 0.00000001490116)
})

