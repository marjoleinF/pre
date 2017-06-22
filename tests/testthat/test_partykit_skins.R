context("Tests the 'skin' versions of the partykit functions")

list.rules <- with(environment(pre), list.rules)

test_that("list.rules gives the same as previously",{
  dat <- airquality[complete.cases(airquality),]
  
  set.seed(2649390)
  tmp <- ctree_setup(Ozone ~ ., data = dat, maxdepth = 3, mtry = Inf)
  mini <- ctree_minimal(tmp$dat, tmp$response, tmp$control, tmp$ytrafo, tmp$terms)
  
  rules <- list.rules(mini)
  
  # save_to_test(rules, "tree_nodes_rules")
  expect_equal(rules, read_to_test("tree_nodes_rules"), tolerance = 0.00000001490116)
})

