context("Tests the 'skin' versions of the partykit functions")

test_that("list.rules gives the same as previously", {
  #####
  # Depth 3 
  dat <- airquality
  
  tree <- ctree(Ozone ~ ., data = dat, 
                control = ctree_control(maxdepth = 3, mtry = Inf))
  
  rules <- list.rules(tree)
  
  expect_length(rules, length(tree) - 2L) # -2 for root and one of the first rules
  
  # save_to_test(rules, "tree_nodes_rules")
  expect_equal(rules, read_to_test("tree_nodes_rules"))
  
  #####
  # Stump
  tree <- ctree(Ozone ~ ., data = dat, 
                control = ctree_control(maxdepth = 1, mtry = Inf))
  
  rules <- list.rules(tree)
  
  expect_length(rules, length(tree) - 2L)
  
  # save_to_test(rules, "tree_nodes_rules_stump")
  expect_equal(rules, read_to_test("tree_nodes_rules_stump"))
})
