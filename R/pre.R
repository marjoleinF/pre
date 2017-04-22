utils::globalVariables("%dopar%")

#' Derive a prediction rule ensemble
#'
#' \code{pre} derives a sparse ensemble of rules and/or linear functions for 
#' prediction of a continuous or binary outcome.
#' 
#' @param formula a symbolic description of the model to be fit of the form 
#' \code{y ~ x1 + x2 + ...+ xn}. If the output variable (left-hand side of the 
#' formala) is a factor, an ensemble for binary classification is created.
#' Otherwise, an ensemble for prediction of a continuous variable is created. 
#' Note that input variables may not have 'rule' as (part of) their name, and 
#' the formula may not exclude the intercept (that is \code{+ 0} or \code{- 1} 
#' may not be used in the right-hand side of the formula).
#' @param data matrix or data.frame containing the variables in the model. When a
#' matrix is specified, it must be of class \code{"numeric"} (the input and output 
#' variable must be continuous; the input variables may be 0-1 coded variables). 
#' When a data.frame is specified, the output variable must be of 
#' class \code{"numeric"} and must be a continuous variable; the input variables 
#' must be of class \code{"numeric"} (for continuous input variables), 
#' \code{"logical"} (for binary variables), \code{"factor"} (for nominal input 
#' variables with 2 or more levels), or \code{"ordered" "factor"} (for 
#' ordered input variables).
#' @param type character. Type of base learners to be included in ensemble. 
#' Defaults to "both" (intial ensemble included both rules and linear functions). 
#' Other option may be "rules" (for prediction rules only) or "linear" (for 
#' linear functions only).
#' @param weights an optional vector of observation weights to be used for 
#' deriving the ensemble.
#' @param sampfrac numeric value greater than 0 and smaller than or equal to 1. 
#' Fraction of randomly selected training observations used to produce each 
#' tree. Values smaller than 1 will result in subsamples being drawn without 
#' replacement (i.e., subsampling), value equal to 1 will result in bootstrap 
#' sampling.
#' @param maxdepth numeric. Maximal number of conditions in rules.
#' @param learnrate numeric. Learning rate for sequentially induced trees. If 
#' \code{NULL} (default), the learnrate is set to .01 for regression and to 0 
#' for classification. Setting the learning rate to values > 0 for classification 
#' dramatically increases computation time.
#' @param removeduplicates logical. Remove rules from the ensemble which have 
#' the exact same support in training data?
#' @param removecomplements logical. Remove rules from the ensemble which have
#' the same support in the training data as the inverse of other rules? 
#' @param mtry numeric. Number of randomly selected predictor variables for 
#' creating each split in each tree. Ignored for nominal output variables if
#' \code{learnrate} > 0.
#' @param thres numeric. Threshold for convergence. 
#' @param standardize logical. Standardize rules and linear terms before 
#' estimating the regression model? As this will also standardize dummy coded
#' factors, users are adviced to use the default: \code{standardize = FALSE}.
#' @param winsfrac numeric. Quantiles of data distribution to be used for 
#' winsorizing linear terms. If set to 0, no winsorizing is performed. Note 
#' that ordinal variables are included as linear terms in estimating the
#' regression model, and will also be winsorized.
#' @param normalize logical. Normalize linear variables before estimating the 
#' regression model? Normalizing gives linear terms the same a priori influence 
#' as a typical rule.
#' @param nfolds numeric. Number of folds to be used in performing cross 
#' validation for determining penalty parameter.
#' @param mod.sel.crit character. Model selection criterion to be used for 
#' deriving the final ensemble. The default is \code{"deviance"}, which uses 
#' squared-error for gaussian models (a.k.a. \code{"mse"}) and binomial deviance 
#' for logistic regression. \code{"class"} would give misclassification error, 
#' \code{"auc"} would give area under the ROC curve. Further, \code{"mse"} or 
#' \code{"mae"} (mean squared and mean absolute error) would measure the deviation 
#' from the fitted mean to the binary or continuous response.
#' @param verbose logical. Should information on the initial and final ensemble 
#' be printed to the command line?
#' @param par.init logical. Should parallel foreach be used to generate initial 
#' ensemble? Only used when \verb{learnrate == 0}. Must register parallel 
#' beforehand, such as doMC or others.
#' @param par.final logical. Should parallel foreach be used to perform cross 
#' validation for selecting the final ensemble? Must register parallel beforehand, 
#' such as doMC or others.
#' @param ntrees numeric. Number of trees to generate for the initial ensemble.
#' @param ... Additional arguments to be passed to 
#' \code{\link[glmnet]{cv.glmnet}}.
#' @note The code for deriving rules from the nodes of trees was taken from an 
#' internal function of the \code{partykit} package of Achim Zeileis and Torsten 
#' Hothorn.
#' @return an object of class \code{pre} 
#' @details Inputs can be continuous, ordered or factor variables. Output can be
#' continuous or binary categorical.
#' @examples \donttest{
#' airq.ens <- pre(Ozone ~ ., data = airquality[complete.cases(airquality),], verbose = TRUE)}
#' @import glmnet partykit datasets
#' @export
pre <- function(formula, data, type = "both", weights = rep(1, times = nrow(data)), 
                sampfrac = .5, maxdepth = 3L, learnrate = NULL, 
                removeduplicates = TRUE, mtry = Inf, ntrees = 500,
                removecomplements = TRUE,
                thres = 1e-07, standardize = FALSE, winsfrac = .025, 
                normalize = TRUE, nfolds = 10L, mod.sel.crit = "deviance", 
                verbose = FALSE, par.init = FALSE, par.final = FALSE, ...)   
{ ###################
  ## Preliminaries ##
  ###################
  
  if(par.init | par.final) {
    if (!("foreach" %in% installed.packages()[,1])) {
      warning("Parallel computation requires package foreach, which is not installed. Argument parallel will be set to FALSE. 
              To run in parallel, download and install package foreach from CRAN, and run again.")   
      par.init <- par.final <- FALSE
    }
  }
  if (!is.data.frame(data)) {
    stop("data should be a data frame.")
  }
  if (length(sampfrac) != 1 || sampfrac < 0.01 || sampfrac > 1) {
    stop("Bad value for 'sampfrac'")
  }
  if (length(type) != 1 || (type != "rules" & type != "both" & type != "linear")) {
    stop("Argument type should equal 'both', 'rules' or 'linear'")
  }
  if (length(winsfrac) != 1 || winsfrac < 0 || winsfrac > 0.5) {
    stop("Bad value for 'winsfrac'.")
  }
  if (!is.logical(verbose)) {
    stop("Bad value for 'verbose'.")
  }  
  orig_data <- data
  data <- model.frame(formula, data, na.action = NULL)
  x_names <- attr(attr(data, "terms"), "term.labels")
  y_name <- names(data)[attr(attr(data, "terms"), "response")]
  formula <- formula(data)
  n <- nrow(data)
  if (is.factor(data[,y_name])) {
    classify <- TRUE
    if (is.null(learnrate)) {
      learnrate <- 0
    }
  } else {
    classify <- FALSE
    if (is.null(learnrate)) {
      learnrate <- .01
    }
  }
  if (!(is.numeric(data[,y_name]) | is.factor(data[,y_name]))) {
    stop("Response variable should be continuous (class numeric) or binary (class 
         factor)")
  }
  if (nlevels(data[,y_name]) > 2) {
    stop("No support for multinomial output variables yet.")
  }
  if (any(sapply(data[,x_names], is.character))) {
    stop("Variables specified in formula and data argument are of class character. 
         Please coerce to class 'numeric', 'factor' or 'ordered' 'factor':", 
         x_names[sapply(data[,x_names], is.character)])
  }
  if (classify & learnrate != 0 & !is.infinite(mtry)) {
    warning("Value specified for mtry will not be used when the outcome variable
            is binary and learnrate > 0", immediate. = TRUE)
  }
  if (any(is.na(data))) {
    data <- data[complete.cases(data),]
    n <- nrow(data)
    warning("Some observations have missing values and will be removed. 
            New sample size is ", n, ".\n", immediate. = TRUE)
  }
  if (verbose) {
    if (classify) {
      cat("A rule ensemble for prediction of a categorical output variable will be 
          created.\n")
    } else {
      cat("A rule ensemble for prediction of a continuous output variable will 
          be created.\n")
    }
  }  
  
  #############################
  ## Derive prediction rules ##
  #############################
  
  if (type != "linear") {
    if (learnrate == 0) { # always use ctree()
      if(par.init) {
        rules <- foreach::foreach(i = 1:ntrees, .combine = "c", .packages = "partykit") %dopar% {
          # Take subsample of dataset
          if (sampfrac == 1) { # then bootstrap
            subsample <- sample(1:n, size = n, replace = TRUE)
          } else { # else subsample
            subsample <- sample(1:n, size = round(sampfrac * n), replace = FALSE)
          }
          subsampledata <- data[subsample,]
          # Make sure ctree() can find object specified by weights argument: 
          environment(formula) <- environment()
          # Grow ctree on subsample:
          tree <- ctree(formula, data = subsampledata, weights = weights[subsample], 
                        maxdepth = maxdepth, mtry = mtry)
          # Collect rules from tree:
          unlist(list.rules(tree))
        }
      } else {
        rules <- c() 
        for(i in 1:ntrees) {
          # Take subsample of dataset
          if (sampfrac == 1) { # then bootstrap
            subsample <- sample(1:n, size = n, replace = TRUE)
          } else { # else subsample
            subsample <- sample(1:n, size = round(sampfrac * n), replace = FALSE)
          }
          subsampledata <- data[subsample,]
          # Make sure ctree() can find object specified by weights argument: 
          environment(formula) <- environment()
          # Grow tree on subsample:
          tree <- ctree(formula, data = subsampledata, weights = weights[subsample], 
                        maxdepth = maxdepth, mtry = mtry)
          # Collect rules from tree:
          rules <- append(rules, unlist(list.rules(tree)))
        }
      }
    }
    if (learnrate > 0) {
      rules <- c()
      if (!classify) {
        y_learn <- data[,y_name]
        for(i in 1:ntrees) {
          # Take subsample of dataset
          if (sampfrac == 1) { # then bootstrap
            subsample <- sample(1:n, size = n, replace = TRUE)
          } else { # else subsample
            subsample <- sample(1:n, size = round(sampfrac * n), replace = FALSE)
          }
          subsampledata <- data[subsample,]
          subsampledata[,y_name] <- y_learn[subsample]
          # Make sure ctree() can find object specified by weights argument: 
          environment(formula) <- environment()
          # Grow tree on subsample:
          tree <- ctree(formula, data = subsampledata, weights = weights[subsample], 
                        maxdepth = maxdepth, mtry = mtry)
          # Collect rules from tree:
          rules <- append(rules, unlist(list.rules(tree)))
          # Substract predictions from current y:
          y_learn <- y_learn - learnrate * predict(tree, newdata = data)
        }
      }
      if (classify) {
        data2 <- data.frame(data, offset = 0)
        glmtreeformula <- formula(paste(paste(y_name, " ~ 1 |"), 
                                        paste(x_names, collapse = "+")))
        for(i in 1:ntrees) {
          # Take subsample of dataset:
          if (sampfrac == 1) { # then bootstrap:
            subsample <- sample(1:n, size = n, replace = TRUE)
          } else { # else subsample:
            subsample <- sample(1:n, size = round(sampfrac * n), replace = FALSE)
          }
          subsampledata <- data2[subsample,]
          # Make sure glmtree() can find object specified by weights argument: 
          environment(formula) <- environment()
          # Grow tree on subsample:
          tree <- glmtree(glmtreeformula, data = subsampledata, family = "binomial",
                          weights = weights[subsample], maxdepth = maxdepth + 1,  
                          offset = offset)
          # Collect rules from tree:
          rules <- append(rules, unlist(list.rules(tree)))
          # Update offset:
          data2$offset <- data2$offset + learnrate * predict(
            tree, newdata = data2, type = "link")
        }
      } 
    }
    nrules <- length(rules)
    if (verbose){
      cat("\nA total of", ntrees, "trees and ", nrules, "rules were 
          generated initially.")
    }
    # Keep unique, non-empty rules only:
    rules <- unique(rules[!rules==""])
    if (verbose) {
      cat("\n\nA total of", nrules - length(rules), "rules were empty
            and removed from the initial ensemble.")
    }
    # Create dataframe with 0-1 coded rules:
    if (length(rules) > 0) {
      rulevars <- data.frame(
        rule1 = as.numeric(with(data, eval(parse(text = rules[[1]])))))
      for(i in 2:length(rules)) {
        rulevars[,paste("rule", i, sep="")] <- as.numeric(
          with(data, eval(parse(text = rules[[i]]))))
      }
      if (removeduplicates) {
        # Remove rules with identical support:
        duplicates <- duplicated(t(rulevars))
        duplicates.removed <- data.frame(name = colnames(rulevars)[duplicates],
                                         description = rules[duplicates])
        rulevars <- rulevars[,!duplicates]
        rules <- rules[!duplicates]
        if (verbose) {
          cat("\n\nA total of", sum(duplicates), "generated rules had 
              support identical to earlier rules and were removed from the initial 
              ensemble ($duplicates.removed shows which, if any).")
        }
      } else {
        duplicates.removed <- NULL
      }
      if (removecomplements) { 
        # remove rules with complement support:
        removed_complement_rules <- c()
        # for rule that has support identical to some earlier rule(s):
        for(i in which(duplicated(apply(rulevars, 2, sd)))) {
          # check whether the rule is a complement of any of the earlier unique rules:
          for(j in 1:i) {
            if (all(rulevars[,i] == (1 - rulevars[,j]))) {
              # add it's name to the list of complement rules:
              removed_complement_rules <- c(removed_complement_rules, names(rulevars)[j])
            }
          }
        }
        # rules with their name in removed_complement_rules should be removed from rulevars and rules
        # and also, some message about the number of rules for which this was the case should be printed if verbose
        complements <- !(names(rulevars) %in% removed_complement_rules)
        rulevars <- rulevars[,!complements]
        rules <- rules[!complements]
        if (verbose) {
          cat("\n\nA total of", length(removed_complement_rules), "generated rules had 
              support that was the complement of the support of earlier rules and were removed from the initial 
              ensemble ($removed_complement_rules shows which, if any).")
        }
      } else {
        removed_complent_rules <- NULL
      }
      if (verbose) {
        cat("\n\nAn initial ensemble consisting of", ncol(rulevars), "rules was 
            succesfully created.")  
      }
    } else {
      warning("No prediction rules could be derived from dataset.", immediate. = TRUE)
    }
  }
  
  ######################################################
  ## Prepare rules, linear terms and outcome variable ##
  ######################################################
  
  x <- data[,x_names]

  # convert ordered categorical predictor variables to linear terms:
  x[,sapply(x, is.ordered)] <- as.numeric(as.character(x[,sapply(x, is.ordered)]))

  if (type == "rules" & length(rules) > 0) {
    x <- rulevars
    x_scales <- NULL
  } else { # if type is not rules, linear terms should be prepared:
    # Winsorize numeric variables (section 5 of F&P(2008)):
    if (winsfrac > 0) {
      wins_points <- data.frame(varname = names(x), value = NA)
      for(i in 1:ncol(x)) {
        if (is.numeric(x[,i])) { 
          lim <- quantile(x[,i], probs = c(winsfrac, 1 - winsfrac))
          x[x[,i] < lim[1], i] <- lim[1]
          x[x[,i] > lim[2], i] <- lim[2]
          wins_points$value[i] <- paste(lim[1], "<=", names(x)[i], "<=", lim[2])
        }
      }
    } else {
      wins_points <- NULL 
    }
    # normalize numeric variables:
    if (normalize) { 
      # Normalize linear terms (section 5 of F&P08), if there are any:
      if (sum(sapply(x, is.numeric)) > 0) {
        x_scales <- sapply(x[sapply(x, is.numeric)], sd, na.rm = TRUE) / 0.4
        x[,sapply(x, is.numeric)] <- scale(x[,sapply(x, is.numeric)], 
                                           center = FALSE, scale = x_scales)
      } else {
        x_scales <- NULL
      }
    } else {
      x_scales <- NULL
    } 
    # If both rules and linear terms are in ensemble, combine both:
    if (type == "both" & length(rules) > 0) {
      x <- data.frame(x, rulevars)
    }
  }
  modmat_formula <- formula(
    paste(" ~ -1 +", paste(colnames(x), collapse = "+")))
  x <- model.matrix(modmat_formula, data = x)
  y <- data[,y_name]
    
  ##################################################
  ## Perform penalized regression on the ensemble ##
  ##################################################
  
  if (classify) {
    family <- "binomial"
  } else {
    family <- "gaussian"
  }
  
  glmnet.fit <- cv.glmnet(x, y, nfolds = nfolds, standardize = standardize, 
                          type.measure = mod.sel.crit, thres = thres, 
                          weights = weights, family = family, parallel = par.final, 
                          ...)
  
  ####################
  ## Return results ##
  ####################
  
  lmin_ind <- which(glmnet.fit$lambda == glmnet.fit$lambda.min)
  l1se_ind <- which(glmnet.fit$lambda == glmnet.fit$lambda.1se)
  if (verbose) {
    cat("\n\nFinal ensemble with minimum cv error: \n  lambda = ", 
        glmnet.fit$lambda[lmin_ind], "\n  number of terms = ", 
        glmnet.fit$nzero[lmin_ind], "\n  mean cv error (se) = ", 
        glmnet.fit$cvm[lmin_ind], " (", glmnet.fit$cvsd[lmin_ind], ")", 
        "\n\nFinal ensemble with cv error within 1se of minimum: \n  lambda = ", 
        glmnet.fit$lambda[l1se_ind],  "\n  number of terms = ", 
        glmnet.fit$nzero[l1se_ind], "\n  mean cv error (se) = ", 
        glmnet.fit$cvm[l1se_ind], " (", glmnet.fit$cvsd[l1se_ind], ")\n", sep="")
  }
  result <- list(glmnet.fit = glmnet.fit, call = match.call(), weights = weights, 
                 data = data, normalize = normalize, x_scales = x_scales, 
                 type = type, x_names = x_names, y_name = y_name, 
                 modmat = x, modmat_formula = modmat_formula, 
                 wins_points = wins_points,
                 classify = classify, formula = formula, orig_data = orig_data)
  if (type != "linear" & length(rules) > 0) {
    result$removed_complement_rules <- removed_complement_rules
    result$duplicates.removed <- duplicates.removed
    result$rules <- data.frame(rule = names(rulevars), description = rules)
    result$rulevars <- rulevars 
  } else {
    result$removed_complement_rules <- NULL
    result$duplicates.removed <- NULL
    result$rules <- NULL
    result$rulevars <- NULL
  }
  class(result) <- "pre"
  return(result)
}



# Internal function for transforming tree into a set of rules:
# Taken and modified from package partykit, written by Achim Zeileis and 
# Torsten Hothorn
list.rules <- function (x, i = NULL, ...) 
{
  if (is.null(i)) 
    i <- partykit::nodeids(x, terminal = TRUE)
  if (length(i) > 1) {
    ret <- sapply(i, list.rules, x = x)
    names(ret) <- if (is.character(i)) 
      i
    else names(x)[i]
    return(ret)
  }
  if (is.character(i) && !is.null(names(x))) 
    i <- which(names(x) %in% i)
  stopifnot(length(i) == 1 & is.numeric(i))
  stopifnot(i <= length(x) & i >= 1)
  i <- as.integer(i)
  dat <- partykit::data_party(x, i)
  if (!is.null(x$fitted)) {
    findx <- which("(fitted)" == names(dat))[1]
    fit <- dat[, findx:ncol(dat), drop = FALSE]
    dat <- dat[, -(findx:ncol(dat)), drop = FALSE]
    if (ncol(dat) == 0) 
      dat <- x$data
  }
  else {
    fit <- NULL
    dat <- x$data
  }
  rule <- c()
  recFun <- function(node) {
    if (partykit::id_node(node) == i) {
      return(NULL)
    }
    kid <- sapply(partykit::kids_node(node), partykit::id_node)
    whichkid <- max(which(kid <= i))
    split <- partykit::split_node(node)
    ivar <- partykit::varid_split(split)
    svar <- names(dat)[ivar]
    index <- partykit::index_split(split)
    if (is.factor(dat[, svar])) {
      if (is.null(index)) 
        index <- ((1:nlevels(dat[, svar])) > partykit::breaks_split(split)) + 
          1
      slevels <- levels(dat[, svar])[index == whichkid]
      srule <- paste(svar, " %in% c(\"", paste(slevels, 
                                               collapse = "\", \"", sep = ""), "\")", sep = "")
    }
    else {
      if (is.null(index)) {
        index <- 1:length(kid)
      }
      breaks <- cbind(c(-Inf, partykit::breaks_split(split)), c(partykit::breaks_split(split), 
                                                      Inf))
      sbreak <- breaks[index == whichkid, ]
      right <- partykit::right_split(split)
      srule <- c()
      if (is.finite(sbreak[1])) {
        srule <- c(srule, paste(svar, ifelse(right, ">", 
                                             ">="), sbreak[1]))
      }
      if (is.finite(sbreak[2])) { 
        srule <- c(srule, paste(svar, ifelse(right, "<=", 
                                             "<"), sbreak[2]))
      }
      srule <- paste(srule, collapse = " & ")
    }
    rule <<- c(rule, srule)
    return(recFun(node[[whichkid]]))
  }
  node <- recFun(partykit::node_party(x))
  paste(rule, collapse = " & ")
}


#' Print method for objects of class pre
#'
#' \code{print.pre} prints information about the generated prediction rule 
#' ensemble to the command line
#' 
#' @param x An object of class \code{\link{pre}}.
#' @param penalty.par.val character. Information for which final prediction rule
#' ensemble(s) should be printed? The ensemble with penalty parameter criterion 
#' yielding minimum cv error (\code{"lambda.min"}) or penalty parameter 
#' yielding error within 1 standard error of minimum cv error 
#' ("\code{lambda.1se}")? Alternatively, a numeric value may be specified, 
#' corresponding to one of the values of lambda in the sequence used by glmnet,
#' for which estimated cv error can be inspected by running \code{x$glmnet.fit}
#' and \code{plot(x$glmnet.fit)}.
#' @param ... Additional arguments, currently not used.
#' @return Prints information about the generated prediction rule ensembles, 
#' @examples \donttest{
#' set.seed(42)
#' airq.ens <- pre(Ozone ~ ., data=airquality[complete.cases(airquality),])
#' coefs <- print(airq.ens)}
#' @export
#' @method print pre
print.pre <- function(x, penalty.par.val = "lambda.1se", ...) {
  if (penalty.par.val == "lambda.1se") {
    lambda_ind <- which(x$glmnet.fit$lambda == x$glmnet.fit$lambda.1se)
    cat("\nFinal ensemble with cv error within 1se of minimum: \n  lambda = ", 
      x$glmnet.fit$lambda[lambda_ind])
  }
  if (penalty.par.val == "lambda.min") {
    lambda_ind <- which(x$glmnet.fit$lambda == x$glmnet.fit$lambda.min)
    cat("Final ensemble with minimum cv error: \n\n  lambda = ", 
        x$glmnet.fit$lambda[lambda_ind])
  }
  if (is.numeric(penalty.par.val)) {
    lambda_ind <- which(round(x$glmnet.fit$lambda, digits = 3) == round(penalty.par.val, digits = 3))
    cat("Final ensemble with lambda = ", round(penalty.par.val, digits = 3))
  }
  cat("\n  number of terms = ", x$glmnet.fit$nzero[lambda_ind], 
      "\n  mean cv error (se) = ", x$glmnet.fit$cvm[lambda_ind], 
        " (", x$glmnet.fit$cvsd[lambda_ind], ") \n\n", sep = "")
  tmp <- coef(x, penalty.par.val = penalty.par.val)
  return(tmp[tmp$coefficient != 0, ])
}





#' Full k-fold cross validation of a pre
#' 
#' \code{cvpre} performs k-fold cross validation on the dataset used to create 
#' the ensemble, providing an estimate of predictive accuracy on future observations.
#' 
#' @param object An object of class \code{\link{pre}}.
#' @param k integer. The number of cross validation folds to be used.
#' @param verbose logical. Should progress of the cross validation be printed 
#' to the command line?
#' @param pclass numeric. Only used for classification. Cut-off value between 
#' 0 and 1 to be used for classifying to second class. 
#' @param penalty.par.val character. Calculate cross-validated error for ensembles 
#' with penalty parameter criterion giving minimum cv error (\code{"lambda.min"}) 
#' or giving cv error that is within 1 standard error of minimum cv error 
#' ("\code{lambda.1se}")? Alternatively, a numeric value may be specified, 
#' corresponding to one of the values of lambda in the sequence used by glmnet,
#' for which estimated cv error can be inspected by running 
#' \code{object$glmnet.fit} and \code{plot(object$glmnet.fit)}.
#' @param parallel logical. Should parallel foreach be used? Must register parallel 
#' beforehand, such as doMC or others.
#' @return A list with three elements: \code{$cvpreds} (a vector with cross-validated
#' predicted y values), \code{$ss} (a vector indicating the cross-validation subsample 
#' each training observation was assigned to) and \code{$accuracy}. For continuous 
#' outputs, accuracy is a list with elements \code{$MSE} (mean squared error on test 
#' observations), \code{$MAE} (mean absolute error on test observations). For 
#' classification, accuracy is a list with elements 
#' \code{$SEL} (mean squared error on predicted probabilities), \code{$AEL} (mean absolute 
#' error on predicted probabilities), \code{$MCR} (average misclassification error rate) 
#' and \code{$table} (table with proportions of (in)correctly classified observations 
#' per class).
#' @examples \donttest{
#' set.seed(42)
#' airq.ens <- pre(Ozone ~ ., data = airquality[complete.cases(airquality),])
#' airq.cv <- cvpre(airq.ens)}
cvpre <- function(object, k = 10, verbose = FALSE, pclass = .5, 
                  penalty.par.val = "lambda.1se", parallel = FALSE) {
  folds <- sample(rep(1:k, length.out = nrow(object$orig_data)), 
               size = nrow(object$orig_data), replace = FALSE)
  if (parallel) {
    cvpreds_unsorted <- foreach::foreach(i = 1:k, .combine = "rbind") %dopar% {
      cl <- object$call
      cl$verbose <- FALSE
      cl$data <- object$orig_data[folds != i,]
      cvobject <- eval(cl)
      if (object$classify) {
        data.frame(fold = rep(i, times = length(folds) - nrow(cvobject$orig_data)), 
                   preds = predict.pre(cvobject, type = "response", 
                                       newdata = object$orig_data[folds == i,], 
                                       penalty.par.val = penalty.par.val))
      } else {
        data.frame(fold = rep(i, times = length(folds) - nrow(cvobject$orig_data)),
                   preds = predict.pre(cvobject, penalty.par.val = penalty.par.val,
                                       newdata = object$orig_data[folds == i,]))
      }
    }
    cvpreds <- rep(NA, times = nrow(object$orig_data))
    for (i in 1:k) {
      cvpreds[folds == i] <- cvpreds_unsorted[cvpreds_unsorted$fold ==i, "preds"]
    }
  } else {
    if (verbose) {
      cat("Running cross validation in fold ")
    }
    cvpreds <- rep(NA, times = nrow(object$orig_data))
    for (i in 1:k) {
      if (verbose) {
        cat(i, " of ", k, ", ", sep = "")
      }
      cl <- object$call
      cl$verbose <- FALSE
      cl$data <- object$orig_data[folds != i,]
      cvobject <- eval(cl)
      if (object$classify) {
        cvpreds[folds == i] <- predict.pre(
          cvobject, newdata = object$orig_data[folds == i,], type = "response", 
          penalty.par.val = penalty.par.val)
      } else {
          cvpreds[folds == i] <- predict.pre(
            cvobject, newdata = object$orig_data[folds == i,], 
            penalty.par.val = penalty.par.val)
      }
      if (verbose & i == k) {
        cat("done!\n")
      }
    }
  }
  accuracy <- list()
  if (object$classify) {
    accuracy$SEL<- mean((as.numeric(object$data[,object$y_name]) - 1 - cvpreds)^2)
    accuracy$AEL <- mean(abs(as.numeric(object$data[,object$y_name]) - 1 - cvpreds))
    cvpreds_d <- as.numeric(cvpreds > .5)
    accuracy$MCR <- 1 - sum(diag(prop.table(table(cvpreds_d, 
                                                  object$data[,object$y_name]))))
    accuracy$table <- prop.table(table(cvpreds_d, object$data[,object$y_name]))
  }
  else {
    accuracy$MSE <- mean((object$data[,object$y_name] - cvpreds)^2)
    accuracy$MAE <- mean(abs(object$data[,object$y_name] - cvpreds))
  }
  result <- list(cvpreds = cvpreds, fold_indicators = folds, accuracy = accuracy)
  return(result)
}





#' Coefficients for the final prediction rule ensemble
#'
#' \code{coef.pre} returns coefficients for prediction rules and linear terms in 
#' the final ensemble
#' 
#' @param object object of class \code{\link{pre}}
#' @param penalty.par.val character. Penalty parameter criterion to be used for 
#' selecting final model: lambda giving minimum cv error (\code{"lambda.min"}) or 
#' lambda giving cv error that is within 1 standard error of minimum cv error 
#' ("\code{lambda.1se}"). Alternatively, a numeric value may be specified, 
#' corresponding to one of the values of lambda in the sequence used by glmnet,
#' for which estimated cv error can be inspected by running 
#' \code{object$glmnet.fit} and \code{plot(object$glmnet.fit)}.
#' @param ... additional arguments to be passed to \code{\link[glmnet]{coef.glmnet}}.
#' @return returns a dataframe with 3 columns: coefficient, rule (rule or 
#' variable name) and description (\code{NA} for linear terms, conditions for 
#' rules).
#' @examples \donttest{
#' set.seed(42)
#' airq.ens <- pre(Ozone ~ ., data=airquality[complete.cases(airquality),])
#' coefs <- coef(airq.ens)}
#' @export
#' @method coef pre
coef.pre <- function(object, penalty.par.val = "lambda.1se", ...)
{
  coefs <- as(coef.glmnet(object$glmnet.fit, s = penalty.par.val, ...), 
              Class = "matrix")
  # coefficients for normalized variables should be unnormalized: 
  if (object$normalize & !is.null(object$x_scales) & object$type != "rules") {
    coefs[names(object$x_scales),] <- coefs[names(object$x_scales),] /
      object$x_scales
  }
  coefs <- data.frame(coefficient = coefs[,1], rule = rownames(coefs))
  if (object$type != "linear" & !is.null(object$rules)) {
    coefs <- merge(coefs, object$rules, all.x = TRUE)
    coefs$description <- as.character(coefs$description)
  } else {
    coefs <- data.frame(rule = coefs$rule, 
                        description = rep(NA, times = nrow(coefs)), 
                        coefficient = coefs[,1])
  }
  if(!is.null(object$wins_points)) { # include winsorizing points in the 
    # description if they were used in generating the ensemble:
    wp <- object$wins_points[!is.na(object$wins_points$value), ]
    coefs[coefs$rule %in% wp$varname, ][
      order(coefs[coefs$rule %in% wp$varname,]$rule), ]$description <- 
      wp[order(wp$varname), ]$value  
  }
  return(coefs[order(abs(coefs$coefficient), decreasing = TRUE),])
}





#' Predicted values based on final unbiased prediction rule ensemble
#'
#' \code{predict.pre} generates predictions based on the final prediction rule
#' ensemble, for training or new (test) observations
#'
#' @param object object of class \code{\link{pre}}.
#' @param newdata optional dataframe of new (test) observations, including all
#' predictor variables used for deriving the prediction rule ensemble.
#' @param penalty.par.val character. Penalty parameter criterion to be used for
#' selecting final model: lambda giving minimum cv error ("lambda.min") or lambda
#' giving cv error that is within 1 standard error of minimum cv error
#' ("lambda.1se"). Alternatively, a numeric value may be specified, 
#' corresponding to one of the values of lambda in the sequence used by glmnet,
#' for which estimated cv error can be inspected by running 
#' \code{object$glmnet.fit} and \code{plot(pbject$glmnet.fit)}.
#' @param type character string. The type of prediction required; the default
#' \code{type = "link"} is on the scale of the linear predictors. Alternatively,
#' for nominal outputs, \code{type = "response"} gives the fitted probabilities
#' and \code{type = "class"} gives the predicted class membership.
#' @param ... currently not used.
#' @details When newdata is not provided, training data included in the specified
#' object is used.
#' @examples \donttest{
#' set.seed(1)
#' train <- sample(1:sum(complete.cases(airquality)), size = 100)
#' set.seed(42)
#' airq.ens <- pre(Ozone ~ ., data = airquality[complete.cases(airquality),][train,])
#' predict(airq.ens)
#' predict(airq.ens, newdata = airquality[complete.cases(airquality),][-train,])}
#' @import Matrix
#' @export
#' @method predict pre
predict.pre <- function(object, newdata = NULL, type = "link",
                         penalty.par.val = "lambda.1se", ...)
{
  if (is.null(newdata)) {
    newdata <- object$modmat
  } else {
    if (!is.data.frame(newdata)) {
      stop("newdata should be a data frame.")
    }
    newdata <- model.frame(object$call$formula, newdata, na.action = NULL)
    # check if newdata has the same columns as object$orig_data:
    if (!all(names(object$data) %in% c(names(newdata), object$y_name))) {
      stop("newdata does not contain all predictor variables from the ensemble")
    } else {
      # take all input variables:
      newdata <- newdata[,names(newdata) %in% object$x_names]
      # add temporary y variable to create model.frame:
      newdata[,object$y_name] <- object$orig_data[,object$y_name][1]
      newdata <- model.frame(object$formula, newdata)
      # check if all variables have the same levels:
      if (!all(unlist(sapply(object$data, levels)) ==
              unlist(sapply(newdata[names(object$data)], levels)))) {
        stop("At least one variable in newdata has different levels than the
             variables used to create the ensemble")
      }
    }
    coefs <- as(coef.glmnet(object$glmnet.fit, s = penalty.par.val),
                Class = "matrix")
    # if there are rules in the ensemble, they should be evaluated:
    if (object$type != "linear") {
      # get names of rules with nonzero and zero coefficients:
      nonzerorulenames <- names(coefs[coefs!=0,])[grep("rule", names(coefs[coefs!=0,]))]
      zerorulenames <- names(coefs[coefs==0,])[grep("rule", names(coefs[coefs==0,]))]
      if (length(nonzerorulenames) > 0) {
        nonzeroterms <- as.character(
          object$rules$description[object$rules$rule %in% nonzerorulenames])
        newrulevars <- data.frame(r1 = as.numeric(with(newdata, eval(parse(
          text = nonzeroterms[1])))))
        names(newrulevars) <- nonzerorulenames[1]
        if (length(nonzerorulenames) > 1) {
          for(i in 2:length(nonzeroterms)) {
            newrulevars[,nonzerorulenames[i]] <- as.numeric(
              with(newdata, eval(parse(text = nonzeroterms[i]))))
          }
        }
        # set all rules with zero coefficients to 0:
        if (length(zerorulenames) > 0) {
          for(i in zerorulenames) {
            newrulevars[,i] <- 0
          }
        }
      } else { # only check and assess rules with non-zero coefficients
        if (length(zerorulenames) > 0) {
          newrulevars <- data.frame(r1 = rep(0, times = nrow(newdata)))
          names(newrulevars) <- zerorulenames[1]
          for(i in zerorulenames[-1]) {
            newrulevars[,i] <- 0
          }
        }
      }
    }
    # convert ordered categorical variables to numeric variables:
    newdata[,sapply(newdata, is.ordered)] <- as.numeric(as.character(
      newdata[,sapply(newdata, is.ordered)]))

    # linear terms normalized before application of glmnet should also be
    # normalized before applying predict.glmnet:
    if (object$normalize & object$type != "rules") {
      newdata[,names(object$x_scales)] <- scale(
        newdata[,names(object$x_scales)], center = FALSE, scale = object$x_scales)
    }
    if (object$type != "linear") {
      newdata <- data.frame(newdata, newrulevars)
    }
    newdata <- MatrixModels::model.Matrix(object$modmat_formula, data = newdata,
                                          sparse = TRUE)
  }

  # Get predictions:
  preds <- predict.cv.glmnet(object$glmnet.fit, newx = newdata, s = penalty.par.val,
                             type = type)[,1]
  return(preds)
}





#' Create partial dependence plot for a single variable
#'
#' \code{singleplot} creates a partial dependence plot, which shows the effect of
#' a predictor variable on the ensemble's predictions
#'
#' @param object an object of class \code{\link{pre}}
#' @param varname character vector of length one, specifying the variable for
#' which the partial dependence plot should be created.
#' penalty.par.val character. Penalty parameter criterion to be used for
#' selecting final model: lambda giving minimum cv error ("lambda.min") or lambda
#' giving cv error that is within 1 standard error of minimum cv error
#' ("lambda.1se"). Alternatively, a numeric value may be specified, 
#' corresponding to one of the values of lambda in the sequence used by glmnet,
#' for which estimated cv error can be inspected by running 
#' \code{object$glmnet.fit} and \code{plot(object$glmnet.fit)}.
#' @param nvals optional numeric vector of length one. For how many values of x
#' should the partial dependence plot be created?
#' @param type character string. Type of prediction to be plotted on y-axis.
#' \code{type = "response"} gives fitted values for continuous outputs and
#' fitted probabilities for nominal outputs. \code{type = "link"} gives fitted
#' values for continuous outputs and linear predictor values for nominal outputs.
#' @param penalty.par.val character. Penalty parameter criterion to be used for
#' selecting final model: lambda giving minimum cv error (\code{"lambda.min"}) or
#' lambda giving cv error that is within 1 standard error of minimum cv error
#' ("\code{lambda.1se}"). Alternatively, a numeric value may be specified, 
#' corresponding to one of the values of lambda in the sequence used by glmnet,
#' for which estimated cv error can be inspected by running 
#' \code{object$glmnet.fit} and \code{plot(object$glmnet.fit)}.
#' @details By default, a partial dependence plot will be created for each unique
#' observed value of the specified predictor variable. When the number of unique
#' observed values is large, this may take a long time to compute. In that case,
#' specifying the nvals argument can substantially reduce computing time. When the
#' nvals argument is supplied, values for the minimum, maximum, and (nvals - 2)
#' intermediate values of the predictor variable will be plotted. Note that nvals
#' can be specified only for numeric and ordered input variables. If the plot is
#' requested for a nominal input variable, the \code{nvals} argument will be
#' ignored and a warning is printed.
#' @examples \donttest{
#' set.seed(42)
#' airq.ens <- pre(Ozone ~ ., data = airquality[complete.cases(airquality),])
#' singleplot(airq.ens, "Temp")}
#' @export
singleplot <- function(object, varname, penalty.par.val = "lambda.1se",
                       nvals = NULL, type = "response")
{
  # preliminaries:
  if (length(varname) != 1) {
    stop("A partial dependence plot should be requested for 1 variable")
  }
  if (!is.character(varname)) {
    stop("Specified varname should be of mode character")
  }
  if (is.factor(object$orig_data[,varname]) & !is.null(nvals)) {
    warning("Plot is requested for variable of class factor. Value specified for
            nvars will be ignored.", immediate. = TRUE)
    nvals <- NULL
  }

  # Generate expanded dataset:
  if (is.null(nvals)) {
    newx <- unique(object$orig_data[,varname])
  } else {
    newx <- seq(
      min(object$orig_data[,varname]), max(object$orig_data[,varname]), length = nvals)
  }
  exp_dataset <- object$orig_data[rep(row.names(object$orig_data), times = length(newx)),]
  exp_dataset[,varname] <- rep(newx, each = nrow(object$orig_data))

  # get predictions:
  exp_dataset$predy <- predict.pre(object, newdata = exp_dataset, type = type,
                                    penalty.par.val = penalty.par.val)

  # create plot:
  plot(aggregate(
    exp_dataset$predy, by = exp_dataset[varname], data = exp_dataset, FUN = mean),
    type = "l", ylab = "predicted y", xlab = varname, main =
      paste("partial dependence on", varname))
  # To be implemented:
  # qntl = trimming factor for plotting numeric variables. Plots are shown for variable values in the range [quantile (qntl) - quantile(1-qntl)]. (Ignored for categorical variables (factors).)
  # nval = maximum number of abscissa evaluation points for numeric variables. (Ignored for categorical variables (factors).)
  # nav = maximum number of observations used for averaging calculations. (larger values provide higher accuracy with a diminishing return; computation grows linearly with nav)
  # catvals = vector of names for values (levels) of categorical variable (factor). (Ignored for numeric variables or length(vars) > 1)
  # samescale = plot vertical scaling flag .
  # samescale = TRUE / FALSE => do/don't require same vertical scale for all plots.
  # horiz = plot orientation flag for categorical variable barplots
  # horiz = T/F => do/don't plot bars horizontally
  # las = label orientation flag for categorical variable plots (horiz = F, only)
  # las = 1 => horizontal orientation of value (level) names stored in catvals (if present)
  # las = 2 => vertical orientation of value (level) names stored in catvals (if present)
  # cex.names = expansion factor for axis names (bar labels) for categorical variable barplots
  # col = color of barplot for categorical variables
  # denqnt = quantile for data density tick marks along upper plot boundary  for numeric variables ( < 1)
  # denqnt <= 0 => no data density tick marks displayed
}





#' Create partial dependence plot for a pair of predictor variables
#'
#' \code{pairplot} creates a partial dependence plot to assess the effects of a
#' pair of predictor variables on the predictions of the ensemble
#'
#' @param object an object of class \code{\link{pre}}
#' @param varnames character vector of length two. Currently, pairplots can only
#' be requested for non-nominal variables. If varnames specifies the name(s) of
#' variables of class \code{"factor"}, an error will be printed.
#' @param penalty.par.val character. Should model be selected with lambda giving
#' minimum cv error ("lambda.min"), or lambda giving cv error that is within 1
#' standard error of minimum cv error ("lambda.1se")? Alternatively, a numeric 
#' value may be specified, corresponding to one of the values of lambda in the 
#' sequence used by glmnet, for which estimated cv error can be inspected by 
#' running \code{object$glmnet.fit} and \code{plot(object$glmnet.fit)}.
#' @param phi numeric. See \code{persp()} documentation.
#' @param theta numeric. See \code{persp()} documentation.
#' @param col character. Optional color to be used for surface in 3D plot.
#' @param nvals optional numeric vector of length 2. For how many values of
#' x1 and x2 should partial dependence be plotted? If \code{NULL}, all observed
#' values for the two predictor variables specified will be used (see details).
#' @param ticktype character string. If \code{"simple"} draws an arrow parallel
#' to the axes to indicate direction of increase; \code{"detailed"} draws ticks 
#' on the axes as in 2D plots.
#' @param nticks the (approximate) number of tick marks to draw on the axes. Has
#' no effect if \code{ticktype = "simple"}.
#' @param type character string. Type of prediction to be plotted on z-axis.
#' \code{type = "response"} gives fitted values for continuous outputs and
#' fitted probabilities for nominal outputs. \code{type = "link"} gives fitted
#' values for continuous outputs and linear predictor values for nominal outputs.
#' @param ... Additional arguments to be passed to \code{\link[graphics]{persp}}.
#' @details By default, partial dependence will be plotted for each combination
#' of 20 values of the specified predictor variables. When \code{nvals = NULL} is
#' specified a dependence plot will be created for every combination of the unique
#' observed values of the two predictor variables specified. Therefore, using
#' \code{nvals = NULL} will often result in long computation times, and / or
#' memory allocation errors. Also, \code{\link{pre}} ensembles derived
#' from training datasets that are very wide or long may result in long
#' computation times and / or memory allocation errors. In such cases, reducing
#' the values supplied to \code{nvals} will reduce computation time and / or
#' memory allocation errors. When the nvals argument is supplied, values for the
#' minimum, maximum, and nvals - 2 intermediate values of the predictor variable
#' will be plotted. Furthermore, if none of the variables specified appears in
#' the final prediction rule ensemble, an error will occur.
#' @note The \code{pairplot} function uses the akima package to construct
#' interpolated surfaces and  has an ACM license that restricts applications
#' to non-commercial usage, see
#' \url{https://www.acm.org/publications/policies/software-copyright-notice}
#' The \code{pairplot} function prints a note refering to this ACM licence.
#' @examples \donttest{
#' set.seed(42)
#' airq.ens <- pre(Ozone ~ ., data = airquality[complete.cases(airquality),])
#' pairplot(airq.ens, c("Temp", "Wind"))}
#' @export
#' @import graphics
pairplot <- function(object, varnames, penalty.par.val = "lambda.1se", phi = 45,
                     theta = 315, col = "cyan", nvals = c(20, 20), ticktype = "detailed",
                     nticks = max(nvals), type = "response", ...)
{
  # preliminaries:
  if (!("akima" %in% installed.packages()[,1])) {
    stop("Function pairplot requires package akima. Download and install package
         akima from CRAN, and run again.")
  }
  if (length(varnames) != 2) {
    stop("Partial dependence should be requested for 2 variables.")
  }
  if (!is.character(varnames)) {
    stop("Specified varname should be of mode character.")
  }
  if (any(sapply(object$orig_data[,varnames], is.factor))) {
    stop("3D partial dependence plots are currently not supported for factors.")
  }
  # generate expanded dataset:
  if (is.null(nvals)){
    newx1 <- unique(object$orig_data[,varnames[1]])
    newx2 <- unique(object$orig_data[,varnames[2]])
  } else {
    newx1 <- seq(min(object$orig_data[,varnames[1]]), max(object$orig_data[,varnames[1]]),
                 length = nvals[1])
    newx2 <- seq(min(object$orig_data[,varnames[2]]), max(object$orig_data[,varnames[2]]),
                 length = nvals[2])
  }
  nobs1 <- length(newx1)
  nobs2 <- length(newx2)
  nobs <- nobs1*nobs2
  exp_dataset <- object$orig_data[rep(row.names(object$orig_data), times = nobs),]
  exp_dataset[,varnames[1]] <- rep(newx1, each = nrow(object$orig_data)*nobs2)
  exp_dataset[,varnames[2]] <- rep(rep(newx2, each = nrow(object$orig_data)),
                                   times = nobs1)

  # get predictions:
  pred_vals <- predict.pre(object, newdata = exp_dataset, type = type,
                            penalty.par.val = penalty.par.val)

  # create plot:
  if (is.null(nvals)) nvals <- 3
  xyz <- akima::interp(exp_dataset[,varnames[1]], exp_dataset[,varnames[2]],
                       pred_vals, duplicate = "mean")
  persp(xyz, xlab = varnames[1], ylab = varnames[2], zlab = "predicted y",
        phi = phi, theta = theta, col = col, ticktype = ticktype,
        nticks = nticks, ...)
  cat("NOTE: function pairplot uses package 'akima', which has an ACM license.
    See also https://www.acm.org/publications/policies/software-copyright-notice.")
}





#' Calculate importances of base learners (rules and linear terms) and input
#' variables
#'
#' \code{importance} calculates importances for rules, linear terms and input
#' variables in the ensemble, and provides a bar plot of variable importances.
#'
#' @param object an object of class \code{\link{pre}}
#' @param plot logical. Should variable importances be plotted?
#' @param ylab character string. Plotting label for y-axis. Only used when
#' \code{plot = TRUE}.
#' @param main character string. Main title of the plot. Only used when
#' \code{plot = TRUE}.
#' @param global logical. Should global importances be calculated? If FALSE,
#' local importances are calculated, given the quantiles of the predictions F(x)
#' in \code{quantprobs}.
#' @param quantprobs optional numeric vector of length two. Only used when
#' \code{global = FALSE} (in which case specification of this argument is still
#' optional). Probabilities for calculating sample quantiles of the range of F(X),
#' over which local importances are calculated. The default provides variable
#' importances calculated over the 25\% highest values of F(X).
#' @param col character string. Plotting color to be used for bars in barplot.
#' @param round integer. Number of decimal places to round numeric results to.
#' If NA (default), no rounding is performed.
#' @param penalty.par.val character. Should model be selected with lambda giving
#' minimum cv error ("lambda.min"), or lambda giving cv error that is within 1
#' standard error of minimum cv error ("lambda.1se")? Alternatively, a numeric 
#' value may be specified, corresponding to one of the values of lambda in the 
#' sequence used by glmnet, for which estimated cv error can be inspected by 
#' running \code{object$glmnet.fit} and \code{plot(object$glmnet.fit)}.
#' @param ... further arguments to be passed to \code{barplot} (only used
#' when \code{plot = TRUE}).
#' @return A list with two dataframes: $baseimps, giving the importances for
#' baselearners in the ensemble, and $varimps, giving the importances for
#' variables that appear and do not appear in the ensemble.
#' @examples \donttest{
#' set.seed(42)
#' airq.ens <- pre(Ozone ~ ., data=airquality[complete.cases(airquality),])
#' # calculate global importances:
#' importance(airq.ens)
#' # calculate local importances (default: over 25% highest predicted values):
#' importance(airq.ens, global = FALSE)
#' # calculate local importances (custom: over 25% lowest predicted values):
#' importance(airq.ens, global = FALSE, quantprobs = c(0, .25))}
#' @export
importance <- function(object, plot = TRUE, ylab = "Importance",
                       main = "Variable importances", global = TRUE,
                       penalty.par.val = "lambda.1se",
                       quantprobs = c(.75, 1), col = "grey", round = NA, ...)
{
  ## Step 1: Calculate the importances of the base learners:

  # get base learner coefficients:
  coefs <- coef.pre(object, penalty.par.val = penalty.par.val)
  # give factors a description:
  coefs$description[is.na(coefs$description)] <-
    paste(as.character(coefs$rule)[is.na(coefs$description)], " ", sep = "")
  coefs <- coefs[order(coefs$rule),]
  # Get sds for every baselearner:
  if (global) {
    sds <- c(0, apply(object$modmat, 2, sd, na.rm = TRUE))
  } else {
    preds <- predict.pre(object, newdata = object$orig_data, type = "response",
                         penalty.par.val = penalty.par.val)
    local_modmat <- object$modmat[preds >= quantile(preds, probs = quantprobs[1]) &
                               preds <= quantile(preds, probs = quantprobs[2]),]
    if (nrow(local_modmat) < 2) {stop("Requested range contains less than 2
                                observations, importances cannot be calculated")}
    sds <- c(0, apply(local_modmat, 2, sd, na.rm = TRUE))
  }
  names(sds)[1] <- "(Intercept)"
  sds <- sds[order(names(sds))]
  if (all(names(sds) != coefs$rule)) {
    stop("There seems to be a problem with the ordering or size of the
         coefficient and sd vectors. Importances cannot be calculated.")
  }

  # baselearner importance is given by abs(coef*st.dev), see F&P section 6):
  baseimps <- data.frame(coefs, sd = sds, imp = abs(coefs$coefficient)*sds)


  ## Step 2: Calculate variable importances:

  # For factors, importances for each level should be added together.
  # first get indicators for assignments in modmat which are not rules:
  inds <- attr(object$modmat, "assign")[-grep("rule", colnames(object$modmat))]
  # add names in modelframe and modelmatrix to baselearner importances:
  frame.mat.conv <- data.frame(
    modmatname = colnames(object$modmat)[-grep("rule", colnames(object$modmat))],
    modframename = attr(attr(object$data, "terms"), "term.labels")[inds])
  baseimps <- merge(frame.mat.conv, baseimps, by.x = "modmatname", by.y = "rule",
                    all.x = TRUE, all.y = TRUE)
  baseimps <- baseimps[baseimps$coefficient != 0,] # helps?
  baseimps <- baseimps[baseimps$description != "(Intercept) ",] # helps?
  # For rules, calculate the number of terms in each rule:
  baseimps$nterms <- NA
  for(i in 1:nrow(baseimps)) {
    # If there is " & " in rule description, there are at least 2 terms/variables 
    # in the base learner:
    if (grepl(" & ", baseimps$description[i])) {
      baseimps$nterms[i] <- length(gregexpr("&", baseimps$description)[[i]]) + 1
    } else {
      baseimps$nterms[i] <- 1 # if not, the number of terms = 1
    }
  }
  # Calculate variable importances:
  varimps <- data.frame(varname = object$x_names, imp = 0)
  # Get importances for rules:
  for(i in 1:nrow(varimps)) { # for every variable:
    # For every baselearner:
    for(j in 1:nrow(baseimps)) {
      # if the variable name appears in the rule:
      #   (Note: EXACT matches are needed, so 1) there should be a space before 
      #     and after the variable name in the rule and thus 2) there should be 
      #     a space added before the description of the rule)
      if(grepl(paste(" ", varimps$varname[i], " ", sep = ""), paste(" ", baseimps$description[j], sep =""))) {
        # then count the number of times it appears in the rule:
        n_occ <- length(gregexpr(paste(" ", varimps$varname[i], " ", sep = ""),
          paste(" ", baseimps$description[j], sep =""), fixed = TRUE)[[1]])
        # and add it to the importance of the variable:
        varimps$imp[i] <- varimps$imp[i] + (n_occ * baseimps$imp[j] /
                                              baseimps$nterms[j])
      }
    }
  }
  # Get importances for factor variables:
  # if the variable appears several times in modframename, add those
  # importances to the variable's importance:
  for(i in object$x_names) {
    if (sum(i == baseimps$modframename, na.rm = TRUE) > 1) {
      varimps$imp[varimps$varname == i] <- sum(varimps$imp[varimps$varname == i],
                      baseimps$imp[i == baseimps$modframename], na.rm = TRUE)
    }
  }

  ## Step 3: return (and plot) importances:
  baseimps <- baseimps[baseimps$imp != 0,]
  baseimps <- baseimps[order(baseimps$imp, decreasing = TRUE),]
  varimps <- varimps[order(varimps$imp, decreasing = TRUE),]
  varimps <- varimps[varimps$imp != 0,]
  if (plot == TRUE & nrow(varimps) > 0) {
    barplot(height = varimps$imp, names.arg = varimps$varname, ylab = ylab,
            main = main, col = col)
  }
  if (!is.na(round)) {
    varimps[,"imp"] <- round(varimps[,"imp"], digits = round)
    baseimps[,c("imp", "coefficient", "sd")] <- round(
      baseimps[,c("imp", "coefficient", "sd")], digits = round)
  }
  return(list(varimps = varimps, baseimps = data.frame(rule = baseimps$modmatname,
    baseimps[baseimps$description != "(Intercept) ", c("description", "imp", "coefficient", "sd")])))
}





#' Compute boostrapped null interaction models
#'
#' \code{bsnullinteract} generates bootstrapped null interaction models,
#' which can be used to derive a reference distribution of the test statistic
#' calculated with \code{\link{interact}}.
#'
#' @param object object of class \code{\link{pre}}.
#' @param nsamp numeric. Number of bootstrapped null interaction models to be
#' derived.
#' @param penalty.par.val character. Which value of the penalty parameter
#' criterion should be used? The value yielding minimum cv error
#' (\code{"lambda.min"}) or penalty parameter yielding error within 1 standard
#' error of minimum cv error ("\code{lambda.1se}")? Alternatively, a numeric 
#' value may be specified, corresponding to one of the values of lambda in the 
#' sequence used by glmnet, for which estimated cv error can be inspected by 
#' running \code{object$glmnet.fit} and \code{plot(object$glmnet.fit)}.
#' @param parallel logical. Should parallel foreach be used to generate initial
#' ensemble? Must register parallel beforehand, such as doMC or others.
#' @param verbose logical. should progress be printed to the command line?
#' @return A list of length \code{nsamp} with null interaction datasets, to be
#' used as input for \code{\link{interact}}.
#' @examples \donttest{
#' set.seed(42)
#' airq.ens <- pre(Ozone ~ ., data=airquality[complete.cases(airquality),])
#' nullmods <- bsnullinteract(airq.ens)}
#' @details Computationally intensive. Progress info is printed to command line.
#' @export
bsnullinteract <- function(object, nsamp = 10, parallel = FALSE,
                           penalty.par.val = "lambda.1se", verbose = FALSE)
{
  # Preliminaries:
  if(parallel) {
    if (!("foreach" %in% installed.packages()[,1])) {
    warning("Parallel computating of function bsnullinteract() requires package foreach,
          which is currently not installed. Argument parallel will be set to FALSE.
          To run in parallel, download and install package foreach from CRAN, and run again.")
      parallel <- FALSE
    }
  }
  # create call for generating bootstrapped null models:
  bsnullmodcall <- object$call
  bsnullmodcall$maxdepth <- 1
  # create call for model allowing for interactions, grown on bootstrapped
  # datasets without interactions:
  bsintmodcall <- object$call
  bsintmodcall$verbose <- FALSE
  # compute boostrapped null datasets (i.e., datasets with no interactions):
  if (parallel) {
    if (verbose) cat("This may take a while.")
    bs.ens <- foreach::foreach(i = 1:nsamp) %dopar% {
      # step 1: Take bootstrap sample {x_p, y_p}:
      bsdataset <- object$orig_data[sample(1:nrow(object$orig_data),
                                           nrow(object$orig_data), replace = TRUE),]
      # step 2: Build F_A, a null interaction model involving main effects only using {x_p, y_p}:
      bsnullmodcall$data <- bsdataset
      bs.ens.null <- eval(bsnullmodcall)
      # step 3: first part of formula 47 of F&P2008:
      # Calculate predictions F_A(x) for original x, using the null interaction model F_A:
      F_a_of_x <- predict.pre(bs.ens.null, newdata = object$orig_data)
      # step 4: third part of formula 47 of F&P2008:
      # Calculate predictions F_A(x_p):
      F_A_of_x_p <- predict.pre(bs.ens.null, newdata = bsdataset,
                            penalty.par.val = penalty.par.val)
      # step 5: Calculate ytilde of formula 47 of F&P2008:
      ytilde <- F_a_of_x + bsdataset[,object$y_name] - F_A_of_x_p
      # step 6: Build a model using (x,ytilde), using the same procedure as was
      # originally applied to (x,y):
      bsintmodcall$data <- object$orig_data
      bsintmodcall$data[,object$y_name] <- ytilde
      eval(match.call(pre, call = bsintmodcall))
    }
  } else {
    bs.ens <- list()
    if (verbose) cat("This may take a while. Computing null model ")
    for(i in 1:nsamp) {
      if (verbose) {cat(i, "of", nsamp, ", ")}
      # step 1: Take bootstrap sample {x_p, y_p}:
      bsdataset <- object$orig_data[sample(1:nrow(object$orig_data),
                                           nrow(object$orig_data), replace = TRUE),]
      # step 2: Build F_A, a null interaction model involving main effects only using {x_p, y_p}:
      bsnullmodcall$data <- bsdataset
      bs.ens.null <- eval(bsnullmodcall)
      # step 3: first part of formula 47 of F&P2008:
      # Calculate predictions F_A(x) for original x, using the null interaction model F_A:
      F_a_of_x <- predict.pre(bs.ens.null, newdata = object$orig_data)
      # step 4: third part of formula 47 of F&P2008:
      # Calculate predictions F_A(x_p):
      F_A_of_x_p <- predict.pre(bs.ens.null, newdata = bsdataset,
                                penalty.par.val = penalty.par.val)
      # step 5: Calculate ytilde of formula 47 of F&P2008:
      ytilde <- F_a_of_x + bsdataset[,object$y_name] - F_A_of_x_p
      # step 6: Build a model using (x,ytilde), using the same procedure as was
      # originally applied to (x,y):
      bsintmodcall$data <- object$orig_data
      bsintmodcall$data[,object$y_name] <- ytilde
      bs.ens[[i]] <- eval(match.call(pre, call = bsintmodcall))
    }
    if (verbose) cat("done!\n")
  }
  return(bs.ens)
}





# Internal function for calculating H statistic (section 8.1, equation 45):
Hsquaredj <- function(object, varname, k = 10, penalty.par.val = NULL, verbose = FALSE) {
  # Calculate the predicted value F(x) of the full model for each observation:
  preds_x <- predict.pre(object, newdata = object$orig_data, penalty.par.val = penalty.par.val)
  # Calculate the expected value of F_j(x_j), over all observed values x_/j,
  # and the expected value of F_/j(x_/j), over all observed values x_j:
  exp_dataset <- object$orig_data[rep(row.names(object$orig_data),
                                 times = nrow(object$orig_data)),]
  exp_dataset[,varname] <- rep(object$orig_data[,varname], each = nrow(object$orig_data))
  # using predict.pre for a hughe dataset may lead to errors, so split
  # computations up in k parts:
  exp_dataset$ids <- sample(1:k, nrow(exp_dataset), replace = TRUE)
  for(i in 1:k) {
    if (verbose) cat(".")
    exp_dataset[exp_dataset$ids==i, "yhat"] <- predict.pre(
      object, newdata = exp_dataset[exp_dataset$ids==i,],
      penalty.par.val = penalty.par.val)
  }
  # expected value of F_j(x_j), over all observed values x_/j:
  exp_dataset$i_xj <- rep(1:nrow(object$orig_data), each = nrow(object$orig_data))
  preds_xj <- aggregate(yhat ~ i_xj, data = exp_dataset, FUN = mean)$yhat
  # expected value of F_/j(x_/j), over all observed values x_j:
  exp_dataset$i_xnotj <-  rep(1:nrow(object$orig_data), times = nrow(object$orig_data))
  preds_xnotj <- aggregate(yhat ~ i_xnotj, data = exp_dataset, FUN = mean)$yhat
  # H should be calculated based on centered functions:
  preds_x <- scale(preds_x, center = TRUE, scale = FALSE)
  preds_xj <- scale(preds_xj, center = TRUE, scale = FALSE)
  preds_xnotj <- scale(preds_xnotj, center = TRUE, scale = FALSE)
  if (sum(preds_x^2) > 0) {
    return(sum((preds_x - preds_xj - preds_xnotj)^2) / sum(preds_x^2))
  }
  if (sum(preds_x^2) == 0) {
    #warning("The denominator for calculating H squared was equal to zero. It was
    #        set to 1e-10 to allow for calculation of H squared", immediate. = TRUE)
    return(sum((preds_x - preds_xj - preds_xnotj)^2) / 1e-10)
  }
}





#' Calculate interaction statistics for user-specified variables
#'
#' \code{interact} calculates test statistics for assessing the strength of
#' interactions between the input variable(s) specified, and all other input
#' variables.
#'
#' @param object an object of class \code{\link{pre}}.
#' @param varnames character vector. Names of variables for which interaction
#' statistics should be calculated. If \code{NULL}, interaction statistics for
#' all predictor variables with non-zeor coefficients will be calculated (which
#' may take a long time).
#' @param k integer. Calculating interaction test statistics is a computationally
#' intensive, so  calculations are split up in several parts to prevent memory
#' allocation errors. If a memory allocation error still occurs, increase k.
#' @param nullmods object with bootstrapped null interaction models, resulting
#' from application of \code{bsnullinteract}.
#' @param penalty.par.val character. Which value of the penalty parameter
#' criterion should be used? The value yielding minimum cv error
#' (\code{"lambda.min"}) or penalty parameter yielding error within 1 standard
#' error of minimum cv error ("\code{lambda.1se}")? Alternatively, a numeric 
#' value may be specified, corresponding to one of the values of lambda in the 
#' sequence used by glmnet, for which estimated cv error can be inspected by 
#' running \code{object$glmnet.fit} and \code{plot(object$glmnet.fit)}.
#' @param parallel logical. Should parallel foreach be used? Must register
#' parallel beforehand, such as doMC or others.
#' @param plot logical Should interaction statistics be plotted?
#' @param col character vector of length two. Color for plotting bars used. Only
#' used when \code{plot = TRUE}. Only first element of vector is used if
#' \code{nullmods = NULL}.
#' @param ylab character string. Label to be used for plotting y-axis.
#' @param main character. Main title for the bar plot.
#' @param legend logical. Should a legend be plotted in the top right corner of the
#' barplot?
#' @param verbose logical. Should progress information be printed to the
#' command line?
#' @param ... Additional arguments to be passed to \code{barplot}.
#' @examples
#' \donttest{
#'  set.seed(42)
#'  airq.ens <- pre(Ozone ~ ., data=airquality[complete.cases(airquality),])
#'  interact(airq.ens, c("Temp", "Wind", "Solar.R"))}
#' @details Can be computationally intensive, especially when nullmods is specified,
#' in which case setting \verb{parallel = TRUE} may improve speed.
#' @return If nullmods is not specified, the function returns the interaction
#' test statistic. If nullmods is specified, the function returns a list,
#' with elements \code{$H}, which is the test statistic of the interaction
#' strength, and \code{$nullH}, which is a vector of test statistics of the
#' interaction in each of the bootstrapped null interaction models. In the barplot,
#' yellow is used for plotting the interaction test statistic. When applicable,
#' blue is used for the mean in the bootstrapped null models.
#' @export
interact <- function(object, varnames = NULL, nullmods = NULL, k = 10, plot = TRUE,
                     penalty.par.val = "lambda.1se", col = c("yellow", "blue"),
                     ylab = "Interaction strength", parallel = FALSE,
                     main = "Interaction test statistics", legend = TRUE,
                     verbose = FALSE, ...)
{ # Preliminaries:
  # Preliminaries:
  if(parallel) {
    if (!("foreach" %in% installed.packages()[,1])) {
      warning("Parallel computating of function bsnullinteract() requires package foreach,
              which is currently not installed. Argument parallel will be set to FALSE.
              To run in parallel, download and install package foreach from CRAN, and run again.")
      parallel <- FALSE
    }
  }
  if (is.null(varnames)) {
    # should only be variables with non-zero importances:
    varnames <- as.character(importance(object, plot = FALSE,
                           penalty.par.val = penalty.par.val)$varimps$varname)
  } else if (!all(varnames %in% object$x_names)) {
    stop("Interaction statistics requested for one or more unknown input variables")
  }
  if (verbose) {
    cat("This will take a while (",
        k * (length(nullmods) + 1) * length(varnames), "dots ). ")
  }
  if (parallel) {
      H <- foreach::foreach(i = 1:length(varnames), .combine = "c") %dopar% {
        # Calculate H_j for the original dataset:
        Hsquaredj(object = object, varname = varnames[i], k = k,
                  penalty.par.val = penalty.par.val, verbose = verbose)
      }
      names(H) <- varnames
      if (!is.null(nullmods)) {
        nullH <- foreach::foreach(i = 1:length(varnames), .combine = "cbind") %dopar% {
          # Calculate H_j for the bootstrapped null models:
          nullH <- c()
          for(j in 1:length(nullmods)) {
           nullH[j] <- Hsquaredj(object = nullmods[[j]], varname = varnames[i],
                                  k = k, penalty.par.val = penalty.par.val,
                                  verbose = verbose)
           }
          nullH
        }
        nullH <- data.frame(nullH)
        names(H) <- colnames(nullH) <- varnames
      }
    } else { # if not parallel computation:
     H <- c()
     if (is.null(nullmods)) {
       for(i in 1:length(varnames)) {
         H[i] <- Hsquaredj(object = object, varname = varnames[i], k = k,
                           penalty.par.val = penalty.par.val, verbose = verbose)
       }
     } else { # Calculate H values for the training data and bootstrapped null models:
       nullH <- data.frame()
       for(i in 1:length(varnames)) {
         H[i] <- Hsquaredj(object = object, varname = varnames[i], k = k,
                           penalty.par.val = penalty.par.val, verbose = verbose)
         for(j in 1:length(nullmods)) {
           nullH[j,i] <- Hsquaredj(object = nullmods[[j]], varname = varnames[i],
                                   k = k, penalty.par.val = penalty.par.val,
                                   verbose = verbose)
         }
       }
       colnames(nullH) <- varnames
    }
    names(H) <- varnames
  }
  if (verbose) cat("\n")
  if (plot) {
    if (is.null(nullmods)) {
      barplot(H, col = col[1], main = main, ...)
    } else {
      nullmeans <- vector()
      for(i in 1:length(varnames)) {
        nullmeans[i] <- mean(nullH[,i])
      }
      H2s <- as.vector(rbind(H, nullmeans))
      barplot(H2s, col = col, ylab = ylab, main = main,
              space = rep_len(1:0, length(H2s)), beside = TRUE,
              names.arg = rep(varnames, each = 2), ...)
       if (legend) {
        legend("topright", c("observed", "bs null mod mean"), bty = "n",
               col = col, pch = 15)
      }
    }
  }
  if(is.null(nullmods)) {
    return(H)
  } else {
    return(list(trainingH2 = H, nullH2 = nullH))
  }
}






#' Plot method for class pre
#'
#' \code{plot.pre} creates one or more plots depicting the rules in the final
#' ensemble as simple decision trees.
#'
#' @param x an object of class \code{\link{pre}}.
#' @param penalty.par.val character. Which value of the penalty parameter
#' criterion should be used? The value yielding minimum cv error
#' (\code{"lambda.min"}) or penalty parameter yielding error within 1 standard
#' error of minimum cv error ("\code{lambda.1se}")? Alternatively, a numeric 
#' value may be specified, corresponding to one of the values of lambda in the 
#' sequence used by glmnet, for which estimated cv error can be inspected by 
#' running \code{x$glmnet.fit} and \code{plot(x$glmnet.fit)}.
#' @param linear.terms logical. Should linear terms be included in the plot?
#' @param nterms numeric. The total number of terms (or rules, if 
#' \code{linear.terms = FALSE}) to be plotted. Default is \code{NULL}, 
#' resulting in all terms of the final ensemble to be plotted.
#' @param max.terms.plot numeric. The maximum number of terms per plot. Rules 
#' are plotted in a square pattern, so \code{is.integer(sqrt(max.terms.plot))} 
#' should return \code{TRUE}, otherwise max.terms.plot will be set to the next 
#' higher value which returns true. The default \code{max.terms.plot = 16} 
#' results in max. 4x4 rules per plot. If the number of terms exceeds the value 
#' specified for max.rules.plot, multiple pages of plots will be created.   
#' @param ask logical. Should user be prompted before starting a new page of
#' plots?
#' @param ... Currently not used.
#' @examples
#' \donttest{
#'  set.seed(42)
#'  airq.ens <- pre(Ozone ~ ., data=airquality[complete.cases(airquality),])
#'  plot(airq.ens)}
#' @export
#' @method plot pre
plot.pre <- function(x, penalty.par.val = "lambda.1se", linear.terms = TRUE, 
                     nterms = NULL, max.terms.plot = 16, ask = FALSE, ...) {
  # Preliminaries:
  if (!("grid" %in% installed.packages()[,1])) {
    stop("Function plot.pre requires package grid. Download and install package
         grid from CRAN, and run again.")
  }
  max.terms.plot <- ceiling(sqrt(max.terms.plot))^2
  # Get nonzero terms from final ensemble:
  nonzeroterms <- importance(x, plot = FALSE, global = TRUE, 
                      penalty.par.val = penalty.par.val)$baseimps
  if (!linear.terms) {
    nonzeroterms <- nonzeroterms[grep("rule", nonzeroterms$rule),]
  }
  if (!is.null(nterms)) {
    nonzeroterms <- nonzeroterms[1:nterms,]
  }
  plot.dim <- rep(min(ceiling(sqrt(nrow(nonzeroterms))), sqrt(max.terms.plot)), times = 2)
  conditions <- list()
  for(i in 1:nrow(nonzeroterms)) { # i is a counter for terms
    if (length(grep("&", nonzeroterms$description[i], )) > 0) { # get rules with multiple conditions:
      conditions[[i]] <- unlist(strsplit(nonzeroterms$description[i], split = " & "))
    } else if (!grepl("rule", nonzeroterms$rule[i])) {
      conditions[[i]] <- "linear" # flag linear terms:
    } else {
      conditions[[i]] <- nonzeroterms$description[i] # get rules with only one condition:
    }
  }
  # Generate a plot for every term:
  for(i in 1:nrow(nonzeroterms)) {
    # track number of plotting page:
    nplot <- floor((i - 1) / max.terms.plot)
    if (conditions[[i]][1] == "linear") { # create plot for linear terms:
      i_plot <- i - nplot * max.terms.plot
      grid::pushViewport(grid::viewport(layout.pos.col = rep(1:plot.dim[2], times = i_plot)[i_plot],
                                        layout.pos.row = ceiling(i_plot/plot.dim[1])))
      grid::grid.text(paste("Linear effect of ", nonzeroterms$rule[i], 
                      "\n\n Importance = ", round(nonzeroterms$imp[i], digits = 3), 
                      sep = ""))
      grid::popViewport()
      if((i-1) %% (plot.dim[1]*plot.dim[2]) == 0) {
        grid::grid.newpage()
        grid::pushViewport(grid::viewport(layout = grid::grid.layout(plot.dim[1], plot.dim[2])))
      }
    } else { # create plot for rules:
      # Create lists of arguments and operators for every condition:
      tmp <- list()
      # check whether ?plot.the operator is " < ", " <= " or "%in%  "
      # split the string using the operator, into the variable name and splitting value, which is used to define split = partysplit(id, value)
      # make it a list:
      for (j in 1:length(conditions[[i]])) {
        condition_j <- conditions[[i]][[j]]
        tmp[[j]] <- character()
        if (length(grep(" > ", condition_j)) > 0) {
          tmp[[j]][1] <- unlist(strsplit(condition_j, " > "))[1]
          tmp[[j]][2] <- " > "
          tmp[[j]][3] <- unlist(strsplit(condition_j, " > "))[2]
        }
        if (length(grep(" <= ", condition_j)) > 0) {
          tmp[[j]][1] <- unlist(strsplit(condition_j, " <= "))[1]
          tmp[[j]][2] <- " <= "
          tmp[[j]][3] <- unlist(strsplit(condition_j, " <= "))[2]
        }
        if (length(grep(" %in% ", condition_j)) > 0) {
          tmp[[j]][1] <- unlist(strsplit(condition_j, " %in% "))[1]
          tmp[[j]][2] <- " %in% "
          tmp[[j]][3] <- unlist(strsplit(condition_j, " %in% "))[2]
        }
      }
      ncond <- length(tmp)
      # generate empty datasets for all the variables appearing in the rules: 
      treeplotdata <- data.frame(matrix(ncol = ncond))
      for (j in 1:ncond) {
        names(treeplotdata)[j] <- tmp[[j]][1]
        if (tmp[[j]][2] == " %in% ") {
          treeplotdata[,j] <- factor(treeplotdata[,j])
          faclevels <- substring(tmp[[j]][3], first = 2)
          faclevels <- gsub(pattern = "\"", replacement = "", x = faclevels, fixed = TRUE)
          faclevels <- gsub(pattern = "(", replacement = "", x = faclevels, fixed = TRUE)
          faclevels <- gsub(pattern = ")", replacement = "", x = faclevels, fixed = TRUE)
          faclevels <- unlist(strsplit(faclevels, ", ",))
          levels(treeplotdata[,j]) <- c(
            levels(x$data[,tmp[[j]][1]])[levels(x$data[,tmp[[j]][1]]) %in% faclevels],
            levels(x$data[,tmp[[j]][1]])[!(levels(x$data[,tmp[[j]][1]]) %in% faclevels)])
          tmp[[j]][3] <- length(faclevels)
        }
      }

      # generate partynode objects for plotting:
      nodes <- list()
      # Construct level 0 of tree (the two terminal nodes), conditional on operator of last condition:
      if (tmp[[ncond]][2] == " > ") { # If condition involves " > ", the tree continues right:
        nodes[[2]] <- list(id = 1L, split = NULL, kids = NULL, surrogates = NULL,
                           info = "exit")
        nodes[[1]] <- list(id = 2L, split = NULL, kids = NULL, surrogates = NULL,
                           info = round(nonzeroterms$coefficient[i], digits = 3))
      } else { # If condition involves " <= " or " %in% " the tree continues left:
        nodes[[2]] <- list(id = 1L, split = NULL, kids = NULL, surrogates = NULL,
                           info = round(nonzeroterms$coefficient[i], digits = 3))
        nodes[[1]] <- list(id = 2L, split = NULL, kids = NULL, surrogates = NULL,
                           info = "exit")
      }
      class(nodes[[1]]) <- class(nodes[[2]]) <- "partynode"

      # if there are > 1 conditions in rule, loop for (nconditions - 1) times:
      if (ncond > 1) {
        for (lev in 1L:(ncond - 1)) { # lev is a counter for the level in the tree
          if (tmp[[lev + 1]][2] == " > ") { # If condition involves " > ", the tree continues right:
            nodes[[lev * 2 + 1]] <- list(id = as.integer(lev * 2 + 1),
                                         split = partysplit(
                                           as.integer(lev), breaks = as.numeric(tmp[[lev]][3])),
                                         kids = list(nodes[[lev * 2 - 1]], nodes[[lev * 2]]),
                                         surrogates = NULL, info = NULL)
            nodes[[lev * 2 + 2]] <- list(id = as.integer(lev * 2 + 2), split = NULL,
                                         kids = NULL, surrogates = NULL, info = "exit")
          } else { # If condition involves " <= " or " %in% " the tree continues left:
            nodes[[lev * 2 + 1]] <- list(id = as.integer(lev * 2 + 1), split = NULL,
                                         kids = NULL, surrogates = NULL, info = "exit")
            nodes[[lev * 2 + 2]] <- list(id = as.integer(lev * 2 + 2),
                                         split = partysplit(
                                           as.integer(lev), breaks = as.numeric(tmp[[lev]][3])),
                                         kids = list(nodes[[lev * 2 - 1]], nodes[[lev * 2]]),
                                         surrogates = NULL, info = NULL)
          }
          class(nodes[[lev * 2 + 1]]) <- class(nodes[[lev * 2 + 2]]) <- "partynode"
        }
      }
      # Construct root node:
      lev <- ncond
      nodes[[lev * 2 + 1]] <- list(id = as.integer(lev * 2 + 2),
                                   split = partysplit(as.integer(lev), breaks = as.numeric(tmp[[lev]][3])),
                                   kids = list(nodes[[lev * 2]], nodes[[lev * 2 - 1]]),
                                   surrogates = NULL, info = NULL)
      class(nodes[[lev * 2 + 1]]) <- "partynode"

      # Plot the rule:
      i_plot <- i - nplot * max.terms.plot
      if((i-1) %% (plot.dim[1]*plot.dim[2]) == 0) {
        grid::grid.newpage()
        grid::pushViewport(grid::viewport(layout = grid::grid.layout(plot.dim[1], plot.dim[2])))
      }
      grid::pushViewport(grid::viewport(layout.pos.col = rep(1:plot.dim[2], times = i_plot)[i_plot],
                                        layout.pos.row = ceiling(i_plot/plot.dim[1])))
      fftree <- party(nodes[[lev * 2 + 1]], data = treeplotdata)
      plot(fftree, newpage = FALSE,
           main = paste(nonzeroterms$rule[i], ": Importance = ", round(nonzeroterms$imp[i], digits = 3), sep = ""),
           inner_panel = node_inner(fftree, id = FALSE),
           terminal_panel = node_terminal(fftree, id = FALSE))
      grid::popViewport()
    }
  }
  if (ask) {
    grDevices::devAskNewPage(ask = FALSE)
  }
}
