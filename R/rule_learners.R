#############################
##
## Rule learner for pre
##
pre_rules <- function(formula, data, weights = rep(1, nrow(data)),
                      y_names, x_names, offset = offset, 
                      learnrate = .01, par.init = FALSE, sampfrac = .5, 
                      mtry = Inf, maxdepth = 3L, ntrees = 500, 
                      tree.control = ctree_control(), use.grad = TRUE, 
                      family = "gaussian", verbose = FALSE, 
                      removeduplicates = TRUE, removecomplements = TRUE,
                      tree.unbiased = TRUE, return.dupl.compl = FALSE, 
                      sparse = FALSE, singleconditions = FALSE, 
                      randomForest = FALSE) {
  
  ## Make sure weights are found by tree-fitting functions.
  ## Most fitting functions search for weights like function lm() does:
  ## weights, subset and offset are evaluated in the same way as variables in formula, 
  ## that is first in data and then in the environment of formula
  environment(formula) <- environment() 
  
  if (randomForest) {
    if (!requireNamespace("randomForest")) {
      stop("Package randomForest has not been installed. Install it and rerun function pre.")
    }
    if (verbose) cat("Fitting randomForest")
    if (mtry > length(x_names)) tree.control$mtry <- length(x_names) ## avoids warning by randomForest
    if (is.null(tree.control$maxnodes)) tree.control$maxnodes <- 2^maxdepth
    rf <- do.call(randomForest::randomForest, 
                  append(list(formula = formula, data = data, weights = weights), tree.control))
    rules <- get_rf_rules(rf, x_names, data, singleconditions = FALSE)
    if (verbose) cat("Done!\n\n")
  } else {
    n <- nrow(data)
    
    ## Prepare glmtree arguments, if necessary:
    if (!use.grad && tree.unbiased) {
      
      ## TODO: Must evaluate tree.control!
      glmtree_args <- tree.control
      glmtree_args$maxdepth <- maxdepth[1] + 1
      glmtree_args$mtry <- mtry
      glmtree_args$formula <- formula(paste(paste(y_names, " ~ 1 |"), 
                                            paste(x_names, collapse = "+")))
      if (!family == "gaussian") {
        glmtree_args$family <- family      
      }
    } else {
      glmtree_args <- NULL
    }
    
    ## Set up subsamples (outside of the loop!):
    if (verbose) cat("\nGenerating random samples...")
    subsample <- if (is.function(sampfrac)) {
      replicate(n = ntrees, sampfrac(n = n, weights = weights))
    } else if (sampfrac == 1) {
      replicate(n = ntrees, sample(1:n, size = n, replace = TRUE, prob = weights))
    } else if (sampfrac < 1) {
      replicate(n = ntrees, sample(1:n, size = round(sampfrac * n), 
                                   replace = FALSE, prob = weights))
    }
    if (is.list(subsample)) {
      ## Add 0s to all samples that do not have maximum length (selecting row 0 return nothing)
      ## and return a
      max_length <- max(sapply(subsample, length))
      subsample <- sapply(subsample, function(x) {
        if (length(x) < max_length) {
          x <- c(x, rep(0L, times = max_length - length(x)))
        }
        x})
    }
    if (verbose) cat(" Done!\n\n")
    
    ## Grow trees
    if (learnrate == 0) {
      
      ## Set up rule learning function:
      fit_tree_return_rules <- function(formula, data, family = NULL, weights,
                                        use.grad = TRUE, tree.unbiased = TRUE,
                                        glmtree_args = NULL, tree.control = NULL) {
        
        ## Make sure weights are found by tree-fitting functions.
        ## Most fitting functions search for weights like e.g. lm() does:
        ## weights, subset and offset are evaluated in the same way as variables in formula, 
        ## that is, first in data and then in the environment of formula
        environment(formula) <- environment()
        
        if (tree.unbiased) {
          if (use.grad) { ## employ ctree
            tree <- partykit::ctree(formula = formula, data = data, 
                                    weights = weights, control = tree.control)
            return(list.rules(tree, removecomplements = removecomplements, 
                              singleconditions = singleconditions))
          } else { ## employ (g)lmtree
            glmtree_args$data <- data
            glmtree_args$weights <- weights          
            if (family == "gaussian") {
              tree <- do.call(partykit::lmtree, args = glmtree_args)
            } else {
              tree <- do.call(partykit::glmtree, args = glmtree_args)
            }
            return(list.rules(tree, removecomplements = removecomplements, 
                              singleconditions = singleconditions))
          }
        } else { # employ rpart
          tree <- rpart::rpart(formula = formula, data = data, 
                               control = tree.control, weights = weights)
          paths <- rpart::path.rpart(tree, nodes = rownames(tree$frame), print.it = FALSE)
          paths <- unname(sapply(sapply(paths, `[`, index = -1), paste, collapse = " & ")[-1])
          if (removecomplements) {
            ## Omit first rule, as it is the complement of a later rule, by definition
            paths <- paths[-1]
          }
          return(paths)
        }
      }
      
      if (par.init) { ## compute in parallel
        rules <- foreach::foreach(i = 1:ntrees, .combine = "c", .packages = c("partykit", "pre")) %dopar% {
          
          if (length(maxdepth) > 1L) {
            if (use.grad) {
              tree.control$maxdepth <- maxdepth[i]
            } else if (tree.unbiased) {
              glmtree_args$maxdepth <- maxdepth[i] + 1L
            }
          }
          fit_tree_return_rules(data = data[subsample[ , i], ], 
                                weights = weights[subsample[ , i]],
                                formula = formula, 
                                family = family, 
                                use.grad = use.grad, 
                                tree.unbiased = tree.unbiased, 
                                glmtree_args = glmtree_args, 
                                tree.control = tree.control)
        }
        
      } else { ## compute in serial
        
        rules <- c()
        if (verbose) {
          cat("Fitting trees\n")
          prog_bar <- txtProgressBar(min = 1, max = ntrees, style = 3)
        }
        for (i in 1:ntrees) {
          
          if (verbose) setTxtProgressBar(prog_bar, i)
          
          if (length(maxdepth) > 1L) {
            if (use.grad) {
              tree.control$maxdepth <- maxdepth[i]
            } else if (tree.unbiased) {
              glmtree_args$maxdepth <- maxdepth[i] + 1L
            }
          }
          rules <- c(rules, 
                     fit_tree_return_rules(data = data[subsample[ , i], ], 
                                           weights = weights[subsample[ , i]],
                                           formula = formula, 
                                           family = family, 
                                           use.grad = use.grad, 
                                           tree.unbiased = tree.unbiased, 
                                           glmtree_args = glmtree_args, 
                                           tree.control = tree.control))
          
        }
      }
      
    } else { ## learnrate > 0
      
      if (verbose) {
        cat("Fitting trees\n")
        prog_bar <- txtProgressBar(min = 1, max = ntrees, style = 3)
      }
      
      rules <- c() ## initialize with empty rule vector
      
      if (use.grad) { ## use ctrees or rpart with y_learn and eta
        
        data_with_y_learn <- data
        ## set initial y and eta value
        if (family == "gaussian") {
          y <- data[[y_names]]
          eta_0 <- weighted.mean(y, weights)
          eta <- rep(eta_0, length(y))
          data_with_y_learn[[y_names]] <- y - eta
        } else if (family == "binomial") {
          y <- data[[y_names]] == levels(data[[y_names]])[1]
          eta_0 <- get_intercept_logistic(y, weights)
          eta <- rep(eta_0, length(y))
          p_0 <- 1 / (1 + exp(-eta))
          data_with_y_learn[[y_names]] <- ifelse(y, log(p_0), log(1 - p_0))
        } else if (family == "poisson") {
          y <- data[[y_names]] 
          eta_0 <- get_intercept_count(y, weights)
          eta <- rep(eta_0, length(y))
          data_with_y_learn[[y_names]] <- y - exp(eta)
        } else if (family == "multinomial") {
          y <- data[y_names]
          ## create dummy variables
          y <- model.matrix(as.formula(paste0(" ~ ", y_names, " - 1")), data = y)
          ## adjust formula used by ctree to involve multiple response variables
          formula <- as.formula(paste(paste(colnames(y), collapse = " + "), "~", 
                                      paste(x_names, collapse = " + ")))
          ## get y_learn:
          eta_0 <- get_intercept_multinomial(y, weights)
          eta <- t(replicate(n = nrow(y), expr = eta_0))
          p_0 <- 1 / (1 + exp(-eta))
          for (i in 1:ncol(y)) {
            y[,i] <- ifelse(y[,i] == 1, log(p_0[,i]), log(1 - p_0[,i]))
          }
          ## omit original response and include dummy-coded response in data
          data_with_y_learn <- cbind(data[ , -which(names(data)== y_names), drop = FALSE], y)
          multinomial_y_names <- names(y)
        } else if (family == "mgaussian") {
          y <- data[ , y_names]
          eta_0 <- apply(y, 2, weighted.mean, weights = rep(1, nrow(y)))
          eta <- t(replicate(n = nrow(y), expr = eta_0))
          data_with_y_learn[,y_names] <- y - eta
        } else if (family == "cox") {
          ## Adjust formula used by ctree and rpart:
          formula <- as.formula(paste0("pseudo_y ~ ", 
                                       paste0(x_names, collapse = " + ")))
          y <- data[ , y_names]
          eta_0 <- 0
          eta <- rep(0, times = nrow(data))
          ngradient_CoxPH <- mboost::CoxPH()@ngradient
          ## omit original response and include pseudo y
          data_with_y_learn <- cbind(data[ , -which(names(data)== y_names), drop = FALSE], y)
          data_with_y_learn$pseudo_y <- ngradient_CoxPH(y = y, f = eta, w = weights)
        }
        
        for(i in 1L:ntrees) {
          
          if (verbose) setTxtProgressBar(prog_bar, i)
          if (length(maxdepth) > 1L) tree.control$maxdepth <- maxdepth[i]
          
          ## Grow tree on subsample
          if (tree.unbiased) {
            tree <- ctree(formula = formula, control = tree.control, 
                          weights = weights[subsample[ , i]],
                          data = data_with_y_learn[subsample[ , i], ])
            ## Collect rules
            rules <- c(rules, list.rules(tree, removecomplements = removecomplements,
                                         singleconditions = singleconditions))
          } else {
            tree <- rpart(formula, control = tree.control, weights = weights[subsample[ , i]],
                          data = data_with_y_learn[subsample[ , i], ])
            paths <- path.rpart(tree, nodes = rownames(tree$frame), print.it = FALSE, pretty = 0)
            paths <- unname(sapply(sapply(paths, `[`, index = -1), paste, collapse = " & ")[-1])
            if (removecomplements) {
              paths <- paths[-1]
            }
            
            rules <- c(rules, paths)
          }
          
          ## Update eta and y_learn
          eta <- eta + learnrate * predict(tree, newdata = data_with_y_learn)
          if (family %in% c("gaussian", "mgaussian")) {
            data_with_y_learn[ , y_names] <- y - eta
          } else if (family == "binomial") {
            data_with_y_learn[[y_names]] <- get_y_learn_logistic(eta, y)
          } else if (family == "poisson") {
            data_with_y_learn[[y_names]] <- get_y_learn_count(eta, y)
          } else if (family == "multinomial") {
            data_with_y_learn[ , multinomial_y_names] <- get_y_learn_multinomial(eta, y)  
          } else if (family == "cox") {
            data_with_y_learn$pseudo_y <- ngradient_CoxPH(y = y, f = eta, w = weights)
          }
        }
        
      } else { ## use.grad is FALSE, employ (g)lmtrees with offset
        
        ## initialize with 0 offset (unless offset, which is a hidden argument, was specified)
        data$.offset <- if (!is.null(offset)) offset else rep(0, times = nrow(data))
        
        for(i in 1:ntrees) {
          
          if (verbose) setTxtProgressBar(prog_bar, i)
          
          ## Take subsample of dataset
          glmtree_args$data <- data[subsample[ , i], ]
          glmtree_args$weights <- weights[subsample[ , i]]
          glmtree_args$offset <- bquote(.offset)
          
          if (length(maxdepth) > 1L) {
            glmtree_args$maxdepth <- maxdepth[i] + 1L
          }
          ## Grow tree on subsample
          if (family == "gaussian") {
            tree <- do.call(lmtree, args = glmtree_args)      
          } else {
            tree <- do.call(glmtree, args = glmtree_args) 
          }
          
          ## Collect rules
          rules <- c(rules, list.rules(tree, removecomplements = removecomplements,
                                       singleconditions = singleconditions))
          ## Update offset (note: do not use a dataset which includes the offset for prediction!!!):
          if (learnrate > 0) {
            if (family == "gaussian") {
              data$.offset <- data$.offset + learnrate * 
                predict(tree, newdata = data)
            } else {
              data$.offset <- data$.offset + learnrate * 
                predict(tree, newdata = data, type = "link")
            }
          }
        }
      }
    }
    if (verbose) try(close(prog_bar), silent = TRUE)
  }
  
  if (length(rules) > 0) {
    ## Keep unique, non-empty rules only
    rules <- unique(rules[!rules==""])
    rules <- rules[!is.na(rules)]
    if (sparse) {
      rules <- .get_most_sparse_rule(rules, data)
    }
    ## Adjust rule format if rpart generated rules
    if (!tree.unbiased && mtry >= length(x_names)) {
      if (any(sapply(data, is.factor))) {
        # replace "=" by " %in% c('"
        for (i in names(data)[sapply(data, is.factor)]) { 
          rules <- gsub(pattern = paste0(i, "="), replacement = paste0(i, " %in% c(\""), 
                       x = rules, fixed = TRUE)
        }
        # replace all "," by "','"
        rules <- gsub(pattern = ",", replacement = "\", \"", x = rules, fixed = TRUE)
        ## add "')" at the end of the string
        rules <- strsplit(x = rules, split = " & ", fixed = TRUE)
        for (i in 1:length(rules)) {
          for (j in names(data)[sapply(data, is.factor)]) {
            if (any(grepl(j, rules[[i]], fixed = TRUE))) {
              rules[[i]][grepl(j, rules[[i]], fixed = TRUE)] <- paste0(
                rules[[i]][grepl(j, rules[[i]], fixed = TRUE)], "\")")
            }        
          }
        }
      }
      rules <- sapply(rules, paste0, collapse = " & ")
      # "<" should be " <" and ">=" should be " >= "
      rules <- gsub(pattern = ">=", replacement = " >= ", fixed = TRUE,
                    x = gsub(pattern = "<", replacement = " <", x = rules, fixed = TRUE))
    }

    if (verbose) {
      cat("\nA total of", ntrees, "trees and ", length(rules), "rules were generated initially.")
    }
  
    rules_obj <- delete_duplicates_complements(
      rules = rules, data = data, 
      removecomplements = removecomplements, 
      removeduplicates = removeduplicates, 
      return.dupl.compl = TRUE, sparse = sparse, 
      keep_rulevars = TRUE)
  
    if (verbose && (removeduplicates || removecomplements)) 
      cat("\n\nA total of", length(rules_obj$duplicates.removed) + length(rules_obj$complements.removed), "generated rules were perfectly collinear with earlier rules and removed from the initial ensemble. \n(fit$duplicates.removed and fit$complements.removed show which, if any).")
    
    if (verbose)
      cat("\n\nAn initial ensemble consisting of", length(rules_obj$rules), "rules was successfully created.")  
    
  } else {
    warn <- "No prediction rules could be derived from dataset. " 
    if (tree.unbiased) {
      warn <- paste0(warn, "Consider increasing the criterion for implementing splits and/or turning off the Bonferroni correction through specification of argument tree.control.")
      if (use.grad) {
        warn <- paste0(warn, "ctree_control(alpha = .5, testtype='Univariate'). ")
      } else {
        warn <- paste0(warn, "mob_control(alpha = .5, bonferroni = FALSE). ")
      }
      warn <- paste0(warn, "(The default for alpha is .05, higher values increase likelihood of splitting.) ")
    }
    warn <- paste0(warn, "Consider increasing the size of samples used for rule generation (by specifying sampfrac=.5 or any other value >.5 and <= 1) in the call to function pre().")
    warning(warn, immediate. = TRUE)
    rules_obj <- list(rules = rules,
                      complements.removed = NULL, 
                      duplicates.removed = NULL, 
                      rulevars = NULL)
  }

  rules_obj
}








#' Get rule learner for gpe which mimics behavior of pre
#'
#' \code{gpe_rules_pre} generates a learner which generates rules like 
#' \code{\link{pre}}, which can be supplied to the \code{\link{gpe}} 
#' base_learner argument.
#' 
#' @inheritParams pre 
#' @examples
#' \donttest{## Obtain same fits with pre and gpe
#' set.seed(42)
#' gpe.mod <- gpe(Ozone ~ ., data = airquality[complete.cases(airquality),],  
#'                base_learners = list(gpe_rules_pre(), gpe_linear()))
#' gpe.mod                
#' set.seed(42)
#' pre.mod <- pre(Ozone ~ ., data = airquality[complete.cases(airquality),],)
#' pre.mod}
#' @export
gpe_rules_pre <- function(learnrate = .01, par.init = FALSE, 
                          mtry = Inf, maxdepth = 3L, ntrees = 500, 
                          tree.control = ctree_control(), use.grad = TRUE, 
                          removeduplicates = TRUE, removecomplements = TRUE,
                          tree.unbiased = TRUE) {
  
  cl <- match.call()
  
  ret <- function(formula, data, weights, sample_func, verbose, family) {
    if (!family %in% c("gaussian", "binomial")) {
      warning("gpe_rules supports only gaussian and binomial family")
    }
    if (any(!complete.cases(data))) {
      warning("data contains missing values'")
    }
    data <- model.frame(Formula::as.Formula(formula), data = data, 
                        na.action = NULL)
    pre_rules_args <- list(
      data = data,
      x_names = attr(attr(data, "terms"), "term.labels"),
      y_names = names(data)[attr(attr(data, "terms"), "response")],
      formula = formula(data), # expands dots in formula
      sampfrac = sample_func,
      weights = if (is.null(cl$weights)) {rep(1L, times = nrow(data))} else {cl$weights},
      verbose = ifelse(is.null(cl$verbose), FALSE, cl$verbose), 
      
      learnrate = ifelse(is.null(cl$learnrate), .01, cl$learnrate), 
      par.init = ifelse(is.null(cl$par.init), FALSE, cl$par.init), 
      mtry = ifelse(is.null(cl$mtry), Inf, cl$mtry), 
      maxdepth = ifelse(is.null(cl$maxdepth), 3L, cl$maxdepth), 
      ntrees = ifelse(is.null(cl$ntrees), 500L, cl$ntrees), 
      tree.control = if (is.null(cl$tree.control)) {
        if (is.null(cl$tree.unbiased)) {
          if (is.null(cl$use.grad)) {
            tree.control <- ctree_control()
            tree.control$maxdepth <- ifelse(is.null(cl$maxdepth), 3L, cl$maxdepth)
            tree.control$mtry <- ifelse(is.null(cl$mtry), Inf, cl$mtry)
            tree.control
          } else if (cl$use.grad) {
            tree.control <- ctree_control()
            tree.control$maxdepth <- ifelse(is.null(cl$maxdepth), 3L, cl$maxdepth)
            tree.control$mtry <- ifelse(is.null(cl$mtry), Inf, cl$mtry)
            tree.control
          } else {
            tree.control <- mob_control()
            tree.control$maxdepth <- 1 + ifelse(is.null(cl$maxdepth), 3L, cl$maxdepth)
            tree.control$mtry <- ifelse(is.null(cl$mtry), Inf, cl$mtry)
            tree.control
          }
        } else if (cl$tree.unbiased) {
          if (is.null(cl$use.grad)) {
            tree.control <- ctree_control()
            tree.control$maxdepth <- ifelse(is.null(cl$maxdepth), 3L, cl$maxdepth)
            tree.control$mtry <- ifelse(is.null(cl$mtry), Inf, cl$mtry)
            tree.control
          } else if (cl$use.grad) {
            tree.control <- ctree_control()
            tree.control$maxdepth <- if (is.null(cl$maxdepth)) {3L} else {cl$maxdepth}
            tree.control$mtry <- ifelse(is.null(cl$mtry), Inf, cl$mtry)
            tree.control
          } else {
            tree.control <- mob_control()            
            tree.control$maxdepth <- 1 + if (is.null(cl$maxdepth)) {3L} else {cl$maxdepth}
            tree.control$mtry <- ifelse(is.null(cl$mtry), Inf, cl$mtry)
            tree.control
          }
        } else {
          tree.control <- rpart.control()
          tree.control$maxdepth <- if (is.null(cl$maxdepth)) {3L} else {cl$maxdepth}
          tree.control
        }
      } else {
        cl$tree.control
      }, 
      use.grad = ifelse(is.null(cl$use.grad), TRUE, cl$use.grad),
      removeduplicates = ifelse(is.null(cl$removeduplicates), TRUE, cl$removeduplicates), 
      removecomplements = ifelse(is.null(cl$removecomplements), TRUE, cl$removecomplements),
      tree.unbiased = ifelse(is.null(cl$tree.unbiased), TRUE, cl$tree.unbiased), 
      return.dupl.compl = FALSE
    )
    rules_obj <- do.call(pre_rules, args = pre_rules_args)
    paste0("rTerm(", rules_obj$rules, ")")
  }
  
  return(ret)
  
}







###################################################
##
## Rule learner based on randomForest::randomForest
##
get_rf_rules <- \(rf, x_names, data, formula, singleconditions = FALSE) {
  
  tree_list <- list()
  for (i in 1:rf$ntree) {
    tree_list[[i]] <- randomForest::getTree(rf, i, labelVar = TRUE)
  }
  
  facs <- names(which(sapply(data[ , x_names], is.factor))) ## TODO; this might not deal with ordered factors well
  fac_levs <- lapply(data[ , facs, drop = FALSE], levels)
  
  ## function that extracts all rules from a single tree
  list_rules_func <- \(tree) {
    
    ## status: -1 is terminal, 1 is not
    terminal_nodes <- which(tree$status == -1L)
    tree$left_or_right_child <- (1:nrow(tree)) %in% tree$`left daughter`
    tree$left_or_right_child <- c("root", ifelse(tree$left_or_right_child, "left", "right")[-1])
    
    conditions <- list()
    for (i in sort(terminal_nodes, decreasing = TRUE)) {
      
      ## i counts terminal node to extract path from
      ## j counts current node (only need to know if left or right)
      ## parent counts parent of current node (need to get the split point)
      
      j <- i
      parent <- c(which(tree$`left daughter` == j), which(tree$`right daughter` == j))
      condition <- character()
      
      while (length(parent) > 0L) {
        
        ## we travel from the terminal node upwards through the tree, so every next condition
        ## should precede the previous one
        if (tree$left_or_right_child[j] == "left") { ## values <= or %in%
          if (tree$`split var`[parent] %in% facs) {
            levs <- fac_levs[[as.character(tree$`split var`[parent])]]
            condition <- c(paste(tree$`split var`[parent], "%in%", 
                                 paste0("c(",
                                        paste0(
                                          levs[as.logical(intToBits(as.integer(tree$`split point`[parent]))[1:length(levs)])],
                                          collapse = ", "), ")")),
                           condition)
          } else {
            condition <- c(paste(tree$`split var`[parent], "<=", tree$`split point`[parent]),
                           condition)
          }
        } else if (tree$left_or_right_child[j] == "right") { ## values > or not %in%
          if (tree$`split var`[parent] %in% facs) {
            levs <- fac_levs[[as.character(tree$`split var`[parent])]]
            condition <- c(paste(tree$`split var`[parent], "%in%", 
                                 paste0("c(",
                                        paste0(
                                          levs[!as.logical(intToBits(as.integer(tree$`split point`[parent]))[1L:length(levs)])],
                                          collapse = ", "), ")")),
                           condition)
          } else {
            condition <- c(paste(tree$`split var`[parent], ">", tree$`split point`[parent]),
                           condition)
          }
        }
        j <- parent
        parent <- c(which(tree$`left daughter` == j), which(tree$`right daughter` == j))
      }
      conditions[[i]] <- if (length(condition) > 2L) {
        if (singleconditions) {
          c(condition,
            sapply(2L:length(condition), \(x) paste(condition[2L:x], collapse = " & ")))         
        } else {
          sapply(2L:length(condition), \(x) paste(condition[2L:x], collapse = " & "))
        }
      } else {
        condition
      }
    }
    conditions
  }
  
  ## return all rules as a character vector
  rules <- unlist(sapply(tree_list, list_rules_func))
  rules <- rules[rules != ""]
  unique(rules)
  ## TODO: process removecomplements (defaults to TRUE) and/or removeduplicates (defaults to TRUE) here to save computation
  ## (it is done after calling get_rf_rules anyway, but checking rule variables for being duplicates or complements is more expensive)
  
}
