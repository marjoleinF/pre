context("Tests the 'skin' versions of the partykit functions")

test_that("list.rules gives the same as previously",{
  #####
  # Depth 3 
  dat <- airquality[complete.cases(airquality),]
  
  tree <- ctree(Ozone ~ ., data = dat, 
                control = ctree_control(maxdepth = 3, mtry = Inf))
  
  rules <- list.rules(tree)
  
  expect_length(rules, length(tree) - 1L)
  
  # save_to_test(rules, "tree_nodes_rules")
  expect_equal(rules, read_to_test("tree_nodes_rules"))
  
  #####
  # Stump
  tree <- ctree(Ozone ~ ., data = dat, 
                control = ctree_control(maxdepth = 1, mtry = Inf))
  
  rules <- list.rules(tree)
  
  expect_length(rules, length(tree) - 1L)
  
  # save_to_test(rules, "tree_nodes_rules_stump")
  expect_equal(rules, read_to_test("tree_nodes_rules_stump"))
})

test_that("list.all_rules_wo_complements gives a subset of list.rules w/ and w/o factors", {
  #####
  # Depth 3
  dat <- airquality[complete.cases(airquality),]
  
  tree <- ctree(Ozone ~ ., data = dat, 
                control = ctree_control(maxdepth = 3, mtry = Inf))
  
  superset <- list.rules(tree)
  subset <- list.all_rules_wo_complements(tree)
  
  expect_true(all(subset %in% superset))
  expect_length(subset, length(partykit::nodeids(tree, terminal = TRUE)) - 1L)
  
  #####
  # Stumps
  tree <- ctree(Ozone ~ ., data = dat, 
                control = ctree_control(maxdepth = 1, mtry = Inf))
  
  superset <- list.rules(tree)
  subset <- list.all_rules_wo_complements(tree)
  
  expect_true(all(subset %in% superset))
  expect_length(subset, length(partykit::nodeids(tree, terminal = TRUE)) - 1L)
  
  #####
  # Depth 3 w/ factors
  dat$Wind_cut <- cut(dat$Wind, breaks = 4)
  
  tree <- ctree(Ozone ~ Temp + Wind_cut, data = dat, 
                control = ctree_control(maxdepth = 3, mtry = Inf))
  
  superset <- list.rules(tree)
  subset <- list.all_rules_wo_complements(tree)
  
  expect_true(all(subset %in% superset))
  expect_length(subset, length(partykit::nodeids(tree, terminal = TRUE)) - 1L)
})
