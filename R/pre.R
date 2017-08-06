utils::globalVariables("%dopar%")

#' Derive a prediction rule ensemble
#'
#' \code{pre} derives a sparse ensemble of rules and/or linear functions for 
#' prediction of a continuous or binary outcome.
#' 
#' @param formula a symbolic description of the model to be fit of the form 
#' \code{y ~ x1 + x2 + ...+ xn}. Response (left-hand side of the formula) 
#' should be of class numeric or of class factor (with two levels). If the 
#' response is a factor, an ensemble for binary classification is created.
#' Otherwise, an ensemble for prediction of a numeric response is created. If
#' the outcome is a non-negative count, this should additionally be specified 
#' by setting\code{family = "poisson"}. Note that input variables may not have 
#' 'rule' as (part of) their name, and the formula may not exclude the intercept 
#' (that is, \code{+ 0} or \code{- 1} may not be used in the right-hand side of 
#' the formula).
#' @param data data.frame containing the variables in the model. Response must
#' be a factor for binary classification, numeric for (count) regression. Input
#' variables must be of class numeric, factor or ordered factor.
#' @param family character. Specification is required only for non-negative 
#' count responses, by specifying \code{family = "poisson"}. Otherwise, 
#' \code{family = "gaussian"} is employed if response specified in formula
#' is numeric and \code{family = "binomial"} is employed if response is a 
#' binary factor. Note that if \code{family = "poisson"} is specified, 
#' \code{\link[partykit]{glmtree}} with an intercept only models in the nodes 
#' will be employed for inducing trees, instead of \code{\link[partykit]{ctree}}. 
#' Although this yields longer computation times, it also yields better 
#' accuracy for count outcomes. 
#' @param use.grad logical. Should binary outcomes use gradient boosting with 
#' regression trees when \code{learnrate > 0}? That is, use 
#' \code{\link[partykit]{ctree}} as in Friedman (2001), without the line search. 
#' By default set to \code{TRUE}, as this yields shorter
#' computation time. If set to \code{FALSE}, \code{\link[partykit]{glmtree}}
#' with intercept only models in the nodes will be employed. This will yield
#' longer computation times (but may increase the likelihood of detecting
#' interactions). 
#' @param weights an optional vector of observation weights to be used for 
#' deriving the ensemble.
#' @param type character. Specifies type of base learners to be included in the 
#' ensemble. Defaults to \code{"both"} (initial ensemble will include both rules 
#' and linear functions). Other option are \code{"rules"} (prediction 
#' rules only) or \code{"linear"} (linear functions only).
#' @param sampfrac numeric. Takes values \eqn{> 0} and \eqn{\leq 1}, representing 
#' the fraction of randomly selected training observations used to produce each 
#' tree. Values \eqn{< 1} will result in sampling without replacement (i.e., 
#' subsampling), a value of 1 will result in sampling with replacement 
#' (i.e., bootstrapping). 
#' @param maxdepth numeric. Maximum number of conditions that can define a rule.
#' @param learnrate numeric. Learning rate or boosting parameter.
#' @param mtry numeric. Number of randomly selected predictor variables for 
#' creating each split in each tree.
#' @param ntrees numeric. Number of trees to generate for the initial ensemble.
#' @param removeduplicates logical. Remove rules from the ensemble which have 
#' the exact same support in training data?
#' @param removecomplements logical. Remove rules from the ensemble which have
#' the same support in the training data as the inverse of other rules? 
#' @param winsfrac numeric. Quantiles of data distribution to be used for 
#' winsorizing linear terms. If set to 0, no winsorizing is performed. Note 
#' that ordinal variables are included as linear terms in estimating the
#' regression model and will also be winsorized.
#' @param normalize logical. Normalize linear variables before estimating the 
#' regression model? Normalizing gives linear terms the same a priori influence 
#' as a typical rule, by dividing the (winsorized) linear term by 2.5 times its 
#' SD.
#' @param standardize logical. Should rules and linear terms be standardized to
#' have SD equal to 1 before estimating the regression model? This will also 
#' standardize the dummified factors, users are advised to use the default 
#' \code{standardize = FALSE}.
#' @param nfolds numeric. Number of cross-validation folds to be used for 
#' selecting the optimal value of the penalty parameter \eqn{\lambda} in selecting
#' the final ensemble.
#' @param verbose logical. Should information on the initial and final ensemble 
#' be printed to the command line?
#' @param par.init logical. Should parallel foreach be used to generate initial 
#' ensemble? Only used when \verb{learnrate == 0} and \code{family != "poisson"}. 
#' Must register parallel beforehand, such as doMC or others.
#' @param par.final logical. Should parallel foreach be used to perform cross 
#' validation for selecting the final ensemble? Must register parallel beforehand, 
#' such as doMC or others.
#' @param tree.control list with control parameters to be passed to the tree 
#' fitting function, see \code{\link[partykit]{ctree_control}}.
#' @param ... Additional arguments to be passed to 
#' \code{\link[glmnet]{cv.glmnet}}.
#' @details Obervations with missing values will be removed prior to analysis.
#' 
#' In rare cases, duplicated variable names may appear in the model.
#' For example, the first variable is a factor named 'V1' and there are also
#' non-factor variables called 'V10' and/or 'V11' and/or 'V12' (etc). Then for 
#' the binary factor V1, dummy contrast variables will be created, called 
#' 'V10', 'V11', 'V12' (etc). As should be clear from this example, this yields 
#' duplicated variable names, which will yield warnings, errors and incorrect 
#' results. Users should prevent this by renaming variables prior to analysis.
#' @note The code for deriving rules from the nodes of trees was taken from an 
#' internal function of the \code{partykit} package of Achim Zeileis and Torsten 
#' Hothorn.
#' @return an object of class \code{pre}, which contains the initial ensemble of 
#' rules and/or linear terms and the final ensembles for a wide range of penalty
#' parameter values. By default, the final ensemble employed by all of the other
#' methods and functions in package \code{pre} is selected using the 'minimum
#' cross validated error plus 1 standard error' criterion. All functions and 
#' methods also take a \code{penalty.parameter.value} argument, which can be
#' used to select a more or less sparse final ensembles. The 
#' \code{penalty.parameter.value} argument takes values \code{"lambda.1se"} 
#' (the default), \code{"lambda.min"}, or a numeric value. Users can assess 
#' the trade of between sparsity and accuracy provided by every possible value 
#' of the penalty parameter (\eqn{\lambda}) by running \code{object$glmnet.fit} 
#' and \code{plot(object$glmnet.fit)}.
#' @details Inputs can be numeric, ordered or factor variables. Reponse can be
#' a numeric, count or binary categorical variable.
#' @examples \donttest{
#' set.seed(42)
#' airq.ens <- pre(Ozone ~ ., data = airquality[complete.cases(airquality),], verbose = TRUE)}
#' @import glmnet partykit datasets
#' @export
#' @seealso \code{\link{print.pre}}, \code{\link{plot.pre}}, 
#' \code{\link{coef.pre}}, \code{\link{importance}}, \code{\link{predict.pre}}, 
#' \code{\link{interact}}, \code{\link{cvpre}} 
#' @references
#' Friedman, J. H. (2001). Greedy function approximation: a gradient boosting machine. \emph{The Annals of Applied Statistics, 29}(5), 1189-1232.
#' 
pre <- function(formula, data, family = c("gaussian", "binomial", "poisson"),
                use.grad = TRUE, weights, type = "both", sampfrac = .5, 
                maxdepth = 3L, learnrate = .01, mtry = Inf, ntrees = 500, 
                removecomplements = TRUE, removeduplicates = TRUE, 
                winsfrac = .025, normalize = TRUE, standardize = FALSE, 
                nfolds = 10L, verbose = FALSE, par.init = FALSE, 
                par.final = FALSE, tree.control, ...) { 
  
  ###################
  ## Preliminaries ##
  ###################
  
  if (missing(weights)) {weights <- rep(1, times = nrow(data))}
  
  if (missing(tree.control)) {
    tree.control <- ctree_control(maxdepth = maxdepth, mtry = mtry)
  } else {
    tree.control$maxdepth <- maxdepth
    tree.control$mtry <- mtry
  }
  
  if (par.final) {
    if (!("foreach" %in% installed.packages()[,1])) {
      warning("Parallel computation requires package foreach, which is not installed. Argument parallel will be set to FALSE. 
              To run in parallel, download and install package foreach from CRAN, and run again.")   
      par.final <- FALSE
    }
  }
  
  if (!is.data.frame(data)) {stop("Data should be a data frame.")}
  
  if (!(is.function(sampfrac))) {
    if (length(sampfrac) != 1 || sampfrac < 0.01 || sampfrac > 1) {
      stop("Bad value for 'sampfrac'")
    }
  }
  
  if (length(type) != 1 || (type != "rules" & type != "both" & type != "linear")) {
    stop("Argument type should equal 'both', 'rules' or 'linear'")
  }
  
  if (length(winsfrac) != 1 || winsfrac < 0 || winsfrac > 0.5) {
    stop("Bad value for 'winsfrac'.")
  }
  if (!is.logical(verbose)) {stop("Bad value for 'verbose'.")}  
  
  ## prepare model frame:
  orig_data <- data
  data <- model.frame(formula, data, na.action = NULL)
  x_names <- attr(attr(data, "terms"), "term.labels")
  y_name <- names(data)[attr(attr(data, "terms"), "response")]
  formula <- formula(data)
  n <- nrow(data)

  ## check and set correct family:
  if (!(is.numeric(data[,y_name]) | is.factor(data[,y_name]))) {
    stop("Response variable should be of class numeric or factor.")
  } else if (family[1] != "poisson") { # if family is gaussian or binomial:
    family <- ifelse(is.numeric(data[,y_name]), "gaussian", "binomial")
    if (verbose) {
      cat(ifelse(family == "gaussian", 
                 "A rule ensemble for prediction of a numeric response will be created.\n",
                 "A rule ensemble for prediction of a binary categorical response will be created.\n"))
    } 
  } else if (verbose) { # if family is poisson and verbose is TRUE, print message:
    cat("A rule ensemble for prediction of a count response will be created.\n")
  }
  
  if (family == "binomial" && (nlevels(data[,y_name]) > 2)) {
    stop("No support for multinomial responses yet.")
  }
  
  if (any(sapply(data[,x_names], is.character))) {
    stop("Variables specified in formula and data argument are of class character. Coerce to class 'numeric', 'factor' or 'ordered' 'factor':", paste(x_names[sapply(data[,x_names], is.character)], sep = ", "))
  }
  
  if (any(is.na(data))) {
    weights <- weights[complete.cases(data)]
    data <- data[complete.cases(data),]
    n <- nrow(data)
    warning("Some observations have missing values and have been removed. New sample size is ", n, ".\n", immediate. = TRUE)
  }
  
  #############################
  ## Derive prediction rules ##
  #############################
  
  if (type != "linear") {
    if (family == "poisson" || 
        (family == "binomial" && !use.grad && learnrate > 0)) {
      glmtreeformula <- formula(paste(paste(y_name, " ~ 1 |"), 
                                      paste(x_names, collapse = "+")))
    }
    
    if (learnrate == 0) {
    ## if learnrate == 0, parallel computation with ctree can be used:
      if (par.init) {
        rules <- foreach::foreach(i = 1:ntrees, .combine = "c", .packages = "partykit") %dopar% {
          # Take subsample of dataset
          if (sampfrac == 1) { # then bootstrap
            subsample <- sample(1:n, size = n, replace = TRUE, prob = weights)
          } else if (sampfrac < 1) { # else subsample
            subsample <- sample(1:n, size = round(sampfrac * n), replace = FALSE, 
                                prob = weights)
          }
          # Grow tree on subsample:
          if (family != "poisson") {
            tree <- ctree(formula = formula, data = data[subsample, ], 
                          control = tree.control)
          } else {
            tree <- glmtree(glmtreeformula, data = data[subsample, ], 
                            family = family, maxdepth = maxdepth + 1, 
                            mtry = mtry)
          }
          # Collect rules from tree:
          rules <- c(rules, list.rules(tree))
        }
      ## if (learnrate == 0 && !par.init), use ctree in a standard for loop:
      } else { # do not compute in parallel:
        rules <- c()
        for(i in 1:ntrees) {
          # Take subsample of dataset
          if (sampfrac == 1) { # then bootstrap
            subsample <- sample(1:n, size = n, replace = TRUE, prob = weights)
          } else if (sampfrac < 1) { # else subsample
            subsample <- sample(1:n, size = round(sampfrac * n), replace = FALSE, 
                                prob = weights)
          }
          # Grow tree on subsample:
          if (family == "poisson") {
            tree <- glmtree(glmtreeformula, data = data[subsample, ], family = family, 
                            maxdepth = maxdepth + 1, mtry = mtry)
          } else {
            tree <- ctree(formula = formula, data = data[subsample, ], 
                          control = tree.control)
          }
          # Collect rules from tree:
          rules <- c(rules, list.rules(tree))
        }
      }
    }
    
    ## If learnrate > 0, induce trees sequentially (no parallel computation):
    if (learnrate > 0) {
      rules <- c()
      
      if (family == "gaussian" || (family == "binomial" && use.grad)) { ## use ctrees with y_learn:
        data_with_y_learn <- data
        if (family == "binomial") {
          y <- data[[y_name]] == levels(data[[y_name]])[1]
          eta_0 <- get_intercept_logistic(y, weights)
          eta <- rep(eta_0, length(y))
          p_0 <- 1 / (1 + exp(-eta))
          data_with_y_learn[[y_name]] <- ifelse(y, log(p_0), log(1 - p_0))
        }
        for(i in 1:ntrees) {
          # Take subsample of dataset:
          if (sampfrac == 1) { # then bootstrap
            subsample <- sample(1:n, size = n, replace = TRUE, prob = weights)
          } else if (sampfrac < 1) { # else subsample
            subsample <- sample(1:n, size = round(sampfrac * n), replace = FALSE, 
                                prob = weights)
          }
          # Grow tree on subsample:
          tree <- ctree(formula = formula, control = tree.control,
                        data = data_with_y_learn[subsample, ])
          # Collect rules from tree:
          rules <- c(rules, list.rules(tree))
          # Substract predictions from current y:
          if (family == "gaussian") {
            data_with_y_learn[[y_name]] <- data_with_y_learn[[y_name]] - 
              learnrate * predict(tree, newdata = data)
          } else { ## family is binomial
            eta <- eta + learnrate * predict(tree, newdata = data)
            data_with_y_learn[[y_name]] <- get_y_learn_logistic(eta, y)
          }
        }
        
      } else { ## use glmtrees with offset: 
        data_with_offset <- data.frame(data, offset = 0)
        for(i in 1:ntrees) {
          # Take subsample of dataset:
          if (sampfrac == 1) { # then bootstrap
            subsample <- sample(1:n, size = n, replace = TRUE, prob = weights)
          } else if (sampfrac < 1) { # else subsample
            subsample <- sample(1:n, size = round(sampfrac * n), replace = FALSE, 
                                prob = weights)
          }
          subsampledata <- data_with_offset[subsample,]
          # Grow tree on subsample:
          tree <- glmtree(glmtreeformula, data = subsampledata, family = family, 
                          maxdepth = maxdepth + 1, mtry = mtry, offset = offset)
          # Collect rules from tree:
          rules <- c(rules, list.rules(tree))
          # Update offset (note that dataset without offset is, and should be, employed for prediction):
          if (learnrate > 0) {
            data_with_offset$offset <- data_with_offset$offset + 
              learnrate * predict(tree, newdata = data, type = "link")
          }
        }
      }
    }
    
    # Keep unique, non-empty rules only:
    rules <- unique(rules[!rules==""])
    if (verbose) {
      cat("\nA total of", ntrees, "trees and ", length(rules), "rules were generated initially.")
    }
    # Create dataframe with 0-1 coded rules:
    if (length(rules) > 0) {
      n_rules <- length(rules)
      rulevars <- matrix(
        NA, nrow = nrow(data), ncol = n_rules, 
        dimnames = list(NULL, paste0("rule", 1:n_rules)))
      names(rules) <- colnames(rulevars)
      
      for(i in 1:n_rules) {
        rulevars[, i] <- with(data, eval(parse(text = rules[[i]])))
      }
      
      if (removeduplicates) {
        # Remove rules with identical support:
        duplicates <- duplicated(rulevars, MARGIN = 2)
        duplicates.removed <- data.frame(name = colnames(rulevars)[duplicates],
                                         description = rules[duplicates])
        rulevars <- rulevars[, !duplicates, drop = FALSE]
        rules <- rules[!duplicates]
      }
      
      if (removecomplements) { 
        # remove rules with complement support:
        sds <- apply(rulevars, 2, sd)
        sds_distinct <- 
          sapply(unique(sds), function(x) c(x, sum(sds == x)))
        
        complements <- vector(mode = "logical", length(sds))
        for(i in 1:ncol(sds_distinct)){
          if(sds_distinct[2, i] < 2)
            next
          
          indices <- which(sds == sds_distinct[1, i])
          for(j in 2:length(indices)){
            indices_prev <- indices[1:(j - 1)] 
            complements[indices_prev] <- 
              complements[indices_prev] | apply(
                rulevars[, indices_prev, drop = F] != rulevars[, indices[j]], 2, all)
          }
        }
        
        complements <- which(complements)
        complements.removed <- data.frame(name = colnames(rulevars)[complements],
                                          description = rules[complements])
        if(length(complements) > 0){
          rulevars <- rulevars[, -complements, drop = FALSE]
          rules <- rules[-complements]
        }
      }
      
      if(!exists("complements.removed"))
        complements.removed <- NULL
      if(!exists("duplicates.removed"))
        duplicates.removed <- NULL
      
      if (verbose && (removeduplicates|| removecomplements)) {
        cat("\n\nA total of", sum(duplicates) + length(complements.removed), "generated rules were perfectly collinear with earlier rules and removed from the initial ensemble. \n($duplicates.removed and $complements.removed show which, if any).")
      }
      
      if (verbose) {
        cat("\n\nAn initial ensemble consisting of", ncol(rulevars), "rules was succesfully created.")  
      }
      storage.mode(rulevars) <- "integer"
      rulevars <- data.frame(rulevars)
    }
    # again check if rules were generated:
    if (length(rules) == 0) {
      warning("No prediction rules could be derived from dataset.", immediate. = TRUE)
      rules <- rulevars <- NULL
    }
  }
  
  if (type == "linear") {rules <- rulevars <- NULL}
  
  ######################################################
  ## Prepare rules, linear terms and outcome variable ##
  ######################################################
  
  if (type == "rules" && length(rules) == 0) {
    warning("No prediction rules could be derived from dataset.")
    return(NULL)
  }
  
  modmat_data <- get_modmat(
    formula = formula, 
    data = data, 
    rules = rules, 
    type = type, 
    winsfrac = winsfrac, 
    x_names = x_names, 
    normalize = normalize)
  y <- modmat_data$y
  x <- modmat_data$x
  x_scales <- modmat_data$x_scales
  modmat_formula <- modmat_data$modmat_formula
  wins_points <- modmat_data$wins_points
  
  # check whether there's duplicates in the variable names:
  # (can happen, for example, due to labeling of dummy indicators for factors)
  if (!(length(unique(colnames(x))) == length(colnames(x)))) { 
    warning("There are variables in the model with overlapping variable names. Rename variables and rerun the analysis. See 'Details' under ?pre.") 
  } 
  
  ##################################################
  ## Perform penalized regression on the ensemble ##
  ##################################################
  
  glmnet.fit <- cv.glmnet(x, y, nfolds = nfolds, weights = weights, 
                          family = family, parallel = par.final, 
                          standardize = standardize, ...)
  
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
                 family = family, formula = formula, orig_data = orig_data)
  if (type != "linear" & length(rules) > 0 & length(names(rulevars)) > 0) {
    result$complements.removed <- complements.removed
    result$duplicates.removed <- duplicates.removed
    result$rules <- data.frame(rule = names(rulevars), description = rules, 
                               stringsAsFactors = FALSE)
    result$rulevars <- rulevars 
  } 
  class(result) <- "pre"
  return(result)
}


get_modmat <- function(
  # Pass these if you already have an object
  modmat_formula = NULL, wins_points = NULL, x_scales = NULL,
  # These should be passed in all calls
  formula, data, rules, type, winsfrac, x_names, normalize){
  if(miss_modmat_formula <- is.null(modmat_formula)) {
    #####
    # Need to define modemat
    str_terms <- if(type != "rules") x_names else character()
    if(type != "linear"){
      str_terms <- c(str_terms, paste0("I(", rules, ")"))
    }
    
    modmat_formula <- paste0(
      ". ~ ", paste0(
        str_terms, collapse = " + "))
    modmat_formula <- update(formula, modmat_formula)
  }
  
  # convert ordered categorical predictor variables to linear terms:
  data[,sapply(data, is.ordered)] <- # Needs to be called on the data.frame
    as.numeric(as.character(data[,sapply(data, is.ordered)]))
  
  data <- model.frame(modmat_formula, data)
  if(miss_modmat_formula)
    modmat_formula <- terms(data) # save terms so model factor levels are keept
  x <- model.matrix(modmat_formula, data = data)
  colnames(x)[
    (ncol(x) - length(rules) + 1):ncol(x)] <- names(rules)
  y <- model.response(data)
  
  #####
  # Remove intercept
  attr_x <- attributes(x)
  attr_x$dimnames[[2]] <- attr_x$dimnames[[2]][-1]
  attr_x$dim[2] <- attr_x$dim[2] - 1
  attr_x$assign <- attr_x$assign[-1]
  x <- x[, colnames(x) != "(Intercept)"]
  
  if(type != "rules"){
    #####
    # if type is not rules, linear terms should be prepared:
    
    # Winsorize numeric variables (section 5 of F&P(2008)):
    if (winsfrac > 0) {
      miss_wins_points <- is.null(wins_points)
      if(miss_wins_points)
        wins_points <- data.frame(varname = x_names, value = NA, lb = NA, ub = NA)
      
      j <- 0
      for(i in x_names) {
        j <- j + 1
        if (is.numeric(data[[i]])) {
          x_idx <- which(
            which(attr(terms(data), "term.labels") == i) == 
              attr_x$assign)
          if(length(x_idx) > 1) # User have made a one to many transformation
            next                # We do not winsorize in this case
          
          if(miss_wins_points){
            lim <- quantile(x[, x_idx], probs = c(winsfrac, 1 - winsfrac))
            wins_points$value[j] <- paste(lim[1], "<=", i, "<=", lim[2])
            wins_points$lb[j] <- lim[1]
            wins_points$ub[j] <- lim[2]
          }
          
          lb <- wins_points$lb[j]
          ub <- wins_points$ub[j]
          
          x[, x_idx][x[, x_idx] < lb] <- lb
          x[, x_idx][x[, x_idx] > ub] <- ub
        }
      }
    }
    
    # normalize numeric variables:
    if (normalize) { 
      # Normalize linear terms (section 5 of F&P08), if there are any:
      needs_scalling <- x_names[sapply(data[x_names], # use data as it is un-transformed 
                                       is.numeric)]
      needs_scalling <- which(colnames(x) %in% x_names)
      if (length(needs_scalling) > 0) {
        if(is.null(x_scales))
          x_scales <- apply(
            x[, needs_scalling, drop = FALSE], 2, sd, na.rm = TRUE) / 0.4
        
        x[, needs_scalling] <- scale(
          x[, needs_scalling, drop = FALSE], center = FALSE, scale = x_scales)
      }
    }
  }
  
  if (!exists("wins_points", inherits = FALSE)) {wins_points <- NULL}
  if (!exists("x_scales", inherits = FALSE)) {x_scales <- NULL}
  
  attributes(x) <- attr_x
  
  list(x = x, y = y, modmat_formula = modmat_formula, 
       x_scales = x_scales, wins_points = wins_points)
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
#' @param digits Number of digits to print
#' @param ... Additional arguments, currently not used
#' @return Prints information about the fitted prediction rule ensemble.
#' @details Note that the cv error is estimated with data that was also used 
#' for learning rules and may be too optimistic. Use cvpre() to obtain an 
#' accurate estimate of future prediction error.
#' @examples \donttest{
#' set.seed(42)
#' airq.ens <- pre(Ozone ~ ., data = airquality[complete.cases(airquality),])
#' print(airq.ens)}
#' @export
#' @method print pre
#' @seealso \code{\link{pre}}, \code{\link{plot.pre}}, 
#' \code{\link{coef.pre}}, \code{\link{importance}}, \code{\link{predict.pre}}, 
#' \code{\link{interact}}, \code{\link{cvpre}} 
print.pre <- function(x, penalty.par.val = "lambda.1se", 
                      digits = getOption("digits"), ...) {
  # function to round values:
  rf <- function(x)
    signif(x, digits)
  
  if (penalty.par.val == "lambda.1se") {
    lambda_ind <- which(x$glmnet.fit$lambda == x$glmnet.fit$lambda.1se)
    cat("\nFinal ensemble with cv error within 1se of minimum: \n  lambda = ", 
        rf(x$glmnet.fit$lambda[lambda_ind]))
  }
  if (penalty.par.val == "lambda.min") {
    lambda_ind <- which(x$glmnet.fit$lambda == x$glmnet.fit$lambda.min)
    cat("Final ensemble with minimum cv error: \n\n  lambda = ", 
        rf(x$glmnet.fit$lambda[lambda_ind]))
  }
  if (is.numeric(penalty.par.val)) {
    lambda_ind <- which(abs(x$glmnet.fit$lambda - penalty.par.val) == min(abs(
      x$glmnet.fit$lambda - penalty.par.val)))
    cat("Final ensemble with lambda = ", rf(x$glmnet.fit$lambda[lambda_ind]))
  }
  cat("\n  number of terms = ", x$glmnet.fit$nzero[lambda_ind], 
      "\n  mean cv error (se) = ", rf(x$glmnet.fit$cvm[lambda_ind]), 
      " (", rf(x$glmnet.fit$cvsd[lambda_ind]), ")", "\n\n  cv error type : ",
      x$glmnet.fit$name, "\n\n", sep = "")
  coefs <- coef(x, penalty.par.val = penalty.par.val)
  coefs <- coefs[coefs$coefficient != 0, ]
  # always put intercept first:
  is_intercept <- 
    if(is.null(coefs$rule))
      rownames(coefs) == "(Intercept)" else coefs$rule == "(Intercept)"
  coefs <- rbind(coefs[is_intercept,], coefs[!is_intercept,])
  
  print(coefs, print.gap = 2, quote = FALSE, row.names = FALSE, digits = digits)
  invisible(coefs)
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
#' @param pclass numeric. Only used for classification. Cut-off value for the 
#' predicted probabilities that should be used to classify observations to the
#' second class. 
#' @param penalty.par.val character. Calculate cross-validated error for ensembles 
#' with penalty parameter criterion giving minimum cv error (\code{"lambda.min"}) 
#' or giving cv error that is within 1 standard error of minimum cv error 
#' ("\code{lambda.1se}")? Alternatively, a numeric value may be specified, 
#' corresponding to one of the values of lambda in the sequence used by glmnet,
#' for which estimated cv error can be inspected by running 
#' \code{object$glmnet.fit} and \code{plot(object$glmnet.fit)}.
#' @param parallel logical. Should parallel foreach be used? Must register parallel 
#' beforehand, such as doMC or others.
#' @return A list with three objects: \code{$cvpreds} (a vector with cross-validated
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
#' @export
#' @seealso \code{\link{pre}}, \code{\link{plot.pre}}, 
#' \code{\link{coef.pre}}, \code{\link{importance}}, \code{\link{predict.pre}}, 
#' \code{\link{interact}}, \code{\link{print.pre}} 
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
      data.frame(fold = rep(i, times = length(folds) - nrow(cvobject$orig_data)), 
                 preds = predict.pre(cvobject, type = "response", 
                                     newdata = object$orig_data[folds == i,], 
                                     penalty.par.val = penalty.par.val))
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
      cvpreds[folds == i] <- predict.pre(
        cvobject, newdata = object$orig_data[folds == i,], type = "response", 
        penalty.par.val = penalty.par.val)
      if (verbose & i == k) {
        cat("done!\n")
      }
    }
  }
  accuracy <- list()
  if (object$family == "binomial") {
    accuracy$SEL<- c(mean((as.numeric(object$data[,object$y_name]) - 1 - cvpreds)^2),
                     sd((as.numeric(object$data[,object$y_name]) - 1 - cvpreds)^2))
    names(accuracy$SEL) <- c("SEL", "se")    
    accuracy$AEL <- c(mean(abs(as.numeric(object$data[,object$y_name]) - 1 - cvpreds)),
                      sd(abs(as.numeric(object$data[,object$y_name]) - 1 - cvpreds)))
    names(accuracy$AEL) <- c("AEL", "se")     
    cvpreds_d <- as.numeric(cvpreds > .5)
    accuracy$MCR <- 1 - sum(diag(prop.table(table(cvpreds_d, 
                                                  object$data[,object$y_name]))))
    accuracy$table <- prop.table(table(cvpreds_d, object$data[,object$y_name]))
  } else {
    accuracy$MSE <- c(mean((object$data[,object$y_name] - cvpreds)^2),
                      sd((object$data[,object$y_name] - cvpreds)^2)/sqrt(length(cvpreds)))
    names(accuracy$MSE) <- c("MSE", "se")
    accuracy$MAE <- c(mean(abs(object$data[,object$y_name] - cvpreds)),
                      sd(abs(object$data[,object$y_name] - cvpreds))/sqrt(length(cvpreds)))
    names(accuracy$MAE) <- c("MAE", "se")
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
#' @details In rare cases, duplucated variable names may appear in the model.
#' For example, when the first variable is named 'V1' and is a factor, and 
#' there is a variable called 'V10' and/or 'V11' and/or 'V12' (etc), which 
#' is/are numeric. For the binary factor V1, dummy contrast variables were 
#' created to fit the model, called 'V10', 'V11', 'V12' (etc). As should be 
#' clear from this example, this yields replicated variable names, which may
#' yield errors or incorrect results. Users should avoid this situation by
#' renaming the variables prior to the analysis.
#' @examples \donttest{
#' set.seed(42)
#' airq.ens <- pre(Ozone ~ ., data = airquality[complete.cases(airquality),])
#' coefs <- coef(airq.ens)}
#' @export
#' @method coef pre
#' @seealso \code{\link{pre}}, \code{\link{plot.pre}}, 
#' \code{\link{cvpre}}, \code{\link{importance}}, \code{\link{predict.pre}}, 
#' \code{\link{interact}}, \code{\link{print.pre}} 
coef.pre <- function(object, penalty.par.val = "lambda.1se", ...)
{
  coefs <- as(coef.glmnet(object$glmnet.fit, s = penalty.par.val, ...), 
              Class = "matrix")
  # coefficients for normalized variables should be unnormalized: 
  if (object$normalize & !is.null(object$x_scales) & object$type != "rules") {
    coefs[names(object$x_scales),] <- coefs[names(object$x_scales),] /
      object$x_scales
  }
  coefs <- data.frame(coefficient = coefs[,1], rule = rownames(coefs), 
                      stringsAsFactors = FALSE)
  # check whether there's duplicates in the variable names:
  # (can happen, for example, due to labeling of dummy indicators for factors)
  if (!(length(unique(coefs$rule)) == length(coefs$rule))) { 
    replicates_in_variable_names <- TRUE
    warning("There are variables in the model with overlapping variable names. This may result in errors, or results may not be valid. See 'Details' under ?coef.pre.") 
  } else {
    replicates_in_variable_names <- FALSE
  }
  if (object$type != "linear" & !is.null(object$rules)) {
    # We set sort to FALSE to get comparable results across platforms
    coefs <- base::merge.data.frame(coefs, object$rules, all.x = TRUE, sort = FALSE)
    coefs$description <- as.character(coefs$description)
  } else {
    coefs <- data.frame(rule = coefs$rule, 
                        description = rep(NA, times = nrow(coefs)), 
                        coefficient = coefs[,1],
                        stringsAsFactors = FALSE)
  }
  # include winsorizing points in the description if they were used in 
  # generating the ensemble (and if there are no duplicate variable names):  
  if (!is.null(object$wins_points) && !replicates_in_variable_names) { 
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
#' \code{object$glmnet.fit} and \code{plot(object$glmnet.fit)}.
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
#' @seealso \code{\link{pre}}, \code{\link{plot.pre}}, 
#' \code{\link{coef.pre}}, \code{\link{importance}}, \code{\link{cvpre}}, 
#' \code{\link{interact}}, \code{\link{print.pre}} 
predict.pre <- function(object, newdata = NULL, type = "link",
                        penalty.par.val = "lambda.1se", ...)
{
  if (is.null(newdata)) {
    newdata <- object$modmat
  } else {
    if (!is.data.frame(newdata)) {
      stop("newdata should be a data frame.")
    }
    
    # Get coefficients
    coefs <- as(coef.glmnet(object$glmnet.fit, s = penalty.par.val),
                Class = "matrix")
    
    # Get model matrix
    winsfrac <- (object$call)$winsfrac
    if(is.null(winsfrac))
      winsfrac <- formals(pre)$winsfrac
    tmp <- get_modmat(
      modmat_formula = object$modmat_formula, 
      wins_points = object$wins_points, 
      x_scales = object$x_scales, 
      formula = object$formula, 
      data = newdata, 
      rules = structure(
        object$rules$description, 
        names = object$rules$rule), 
      type = object$type, 
      winsfrac = winsfrac,
      x_names = object$x_names, 
      normalize = object$normalize)
    
    newdata <- tmp$x
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
#' @param ... Further arguments to be passed to 
#' \code{\link[graphics]{plot.default}}.
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
#' @seealso \code{\link{pre}}, \code{\link{pairplot}}
singleplot <- function(object, varname, penalty.par.val = "lambda.1se",
                       nvals = NULL, type = "response", ...)
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
    type = "l", ylab = "predicted y", xlab = varname, ...)
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
#' @param type character string. Type of plot to be generated. 
#' \code{type = "heatmap"} yields a heatmap plot, \code{type = "contour"} yields 
#' a contour plot, \code{type = "both"} yields a heatmap plot with added contours,
#' \code{type = "perspective"} yields a three dimensional plot.
#' @param nvals optional numeric vector of length 2. For how many values of
#' x1 and x2 should partial dependence be plotted? If \code{NULL}, all observed
#' values for the two predictor variables specified will be used (see details).
#' @param pred.type character string. Type of prediction to be plotted on z-axis.
#' \code{pred.type = "response"} gives fitted values for continuous outputs and
#' fitted probabilities for nominal outputs. \code{pred.type = "link"} gives fitted
#' values for continuous outputs and linear predictor values for nominal outputs.
#' @param ... Additional arguments to be passed to \code{\link[graphics]{image}}, 
#' \code{\link[graphics]{contour}} or \code{\link[graphics]{persp}} (depending on
#' whether \code{type} is specified to be \code{"heatmap"}, \code{"contour"}, \code{"both"} 
#' or \code{"perspective"}).
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
#' 
#' @note Function \code{pairplot} uses package akima to construct interpolated 
#' surfaces and  has an ACM license that restricts applications to non-commercial 
#' usage, see 
#' \url{https://www.acm.org/publications/policies/software-copyright-notice}
#' Function \code{pairplot} prints a note referring to this ACM licence.
#' @examples \donttest{
#' set.seed(42)
#' airq.ens <- pre(Ozone ~ ., data = airquality[complete.cases(airquality),])
#' pairplot(airq.ens, c("Temp", "Wind"))}
#' @export
#' @import graphics
#' @export
#' @seealso \code{\link{pre}}, \code{\link{singleplot}} 
pairplot <- function(object, varnames, type = "both",
                     penalty.par.val = "lambda.1se", 
                     nvals = c(20, 20), pred.type = "response", ...)
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
  pred_vals <- predict.pre(object, newdata = exp_dataset, type = pred.type,
                           penalty.par.val = penalty.par.val)
  
  # create plot:
  if (is.null(nvals)) {nvals <- 3}
  xyz <- akima::interp(exp_dataset[,varnames[1]], exp_dataset[,varnames[2]],
                       pred_vals, duplicate = "mean")
  if (type == "heatmap" || type == "both") {
    if (is.null(match.call()$col)) {
      image(xyz, xlab = varnames[1], ylab = varnames[2], 
            col = rev(grDevices::heat.colors(12)), ...)
    } else {
      image(xyz, xlab = varnames[1], ylab = varnames[2], ...)
    }
    if (type == "both") {
      contour(xyz, add = TRUE)
    }
  }
  if (type == "contour") {
    contour(xyz, xlab = varnames[1], ylab = varnames[2], ...) 
  }
  if (type == "perspective") {
    persp(xyz, xlab = varnames[1], ylab = varnames[2], zlab = "predicted y", ...)
  }
  message("NOTE: function pairplot uses package 'akima', which has an ACM license. See also https://www.acm.org/publications/policies/software-copyright-notice.")
}


#' Calculate importances of baselearners (rules and linear terms) and input
#' variables
#'
#' \code{importance} calculates importances for rules, linear terms and input
#' variables in the ensemble, and provides a bar plot of variable importances.
#'
#' @param object an object of class \code{\link{pre}}
#' @param standardize logical. Should baselearner importances be standardized 
#' with respect to the outcome variable? If \code{TRUE}, baselearner importances 
#' have a minimum of 0 and a maximum of 1. Only used for ensembles with 
#' numeric (non-count) response variables.
#' @param global logical. Should global importances be calculated? If 
#' \code{FALSE}, local importances will be calculated, given the quantiles 
#' of the predictions F(x) in \code{quantprobs}.
#' @param quantprobs optional numeric vector of length two. Only used when
#' \code{global = FALSE}. Probabilities for calculating sample quantiles of the 
#' range of F(X), over which local importances are calculated. The default 
#' provides variable importances calculated over the 25\% highest values of F(X).
#' @param penalty.par.val character. Should model be selected with lambda yielding
#' minimum cv error ("lambda.min"), or lambda giving cv error that is within 1
#' standard error of minimum cv error ("lambda.1se")? Alternatively, a numeric 
#' value may be specified, corresponding to one of the values of lambda in the 
#' sequence used by glmnet.
#' @param round integer. Number of decimal places to round numeric results to.
#' If NA (default), no rounding is performed.
#' @param plot logical. Should variable importances be plotted?
#' @param ylab character string. Plotting label for y-axis. Only used when
#' \code{plot = TRUE}.
#' @param main character string. Main title of the plot. Only used when
#' \code{plot = TRUE}.
#' @param diag.xlab logical. Should variable names be printed diagonally (that
#' is, in a 45 degree angle)? Alternatively, variable names may be printed 
#' vertically by specifying \code{diag.xlab = FALSE, las = 2}.
#' @param diag.xlab.hor numeric. Horizontal adjustment for lining up variable
#' names with bars in the plot if variable names are printed diagonally.
#' @param diag.xlab.vert positive integer. Vertical adjustment for position
#' of variable names, if printed diagonally. Corresponds to the number of 
#' character spaces added after variable names. 
#' @param cex.axis numeric. The magnification to be used for axis annotation
#' relative to the current setting of \code{cex}.
#' @param ... further arguments to be passed to \code{barplot} (only used
#' when \code{plot = TRUE}).
#' @return A list with two dataframes: \code{$baseimps}, giving the importances 
#' for baselearners in the ensemble, and \code{$varimps}, giving the importances 
#' for all predictor variables.
#' @examples \donttest{
#' set.seed(42)
#' airq.ens <- pre(Ozone ~ ., data = airquality[complete.cases(airquality),])
#' # calculate global importances:
#' importance(airq.ens)
#' # calculate local importances (default: over 25% highest predicted values):
#' importance(airq.ens, global = FALSE)
#' # calculate local importances (custom: over 25% lowest predicted values):
#' importance(airq.ens, global = FALSE, quantprobs = c(0, .25))}
#' @export
#' #' @seealso \code{\link{pre}}
importance <- function(object, standardize = FALSE, global = TRUE,
                       quantprobs = c(.75, 1), penalty.par.val = "lambda.1se", 
                       round = NA, plot = TRUE, ylab = "Importance",
                       main = "Variable importances", diag.xlab = TRUE, 
                       diag.xlab.hor = 0, diag.xlab.vert = 2,
                       cex.axis = 1, ...)
{
  ## Step 1: Calculate the importances of the base learners:
  
  # get base learner coefficients:
  coefs <- coef.pre(object, penalty.par.val = penalty.par.val)
  # only continue when there are nonzero terms besides intercept:
  if (sum(coefs$coefficient != 0) > 1) { 
    # give factors a description:
    coefs$description[is.na(coefs$description)] <-
      paste0(as.character(coefs$rule)[is.na(coefs$description)], " ")
    coefs <- coefs[order(coefs$rule),]
    # Get sds for every baselearner:
    if (global) {
      # object$x_scales should be used to get correct SDs for linear terms:
      sds <- c(0, apply(object$modmat, 2, sd, na.rm = TRUE))  
      if (standardize) {
        sd_y <- sd(object$data[,object$y_name])
      }
      if(object$normalize) {
        sds[names(object$x_scales)] <- sds[names(object$x_scales)] * object$x_scales
      }
    } else {
      preds <- predict.pre(object, newdata = object$orig_data, type = "response",
                           penalty.par.val = penalty.par.val)
      local_modmat <- object$modmat[preds >= quantile(preds, probs = quantprobs[1]) &
                                      preds <= quantile(preds, probs = quantprobs[2]),]
      if (nrow(local_modmat) < 2) {stop("Selected subregion contains less than 2
                                        observations, importances cannot be calculated")}
      # object$x_scales should be used to get correct SDs for linear terms:
      sds <- c(0, apply(local_modmat, 2, sd, na.rm = TRUE))
      if(object$normalize) {
        sds[names(object$x_scales)] <- sds[names(object$x_scales)] * object$x_scales
      }
      if (standardize) {
        sd_y <- sd(object$data[preds >= quantile(preds, probs = quantprobs[1]) & 
                                 preds <= quantile(preds, probs = quantprobs[2]),
                               object$y_name])
      }
    }
    names(sds)[1] <- "(Intercept)"
    sds <- sds[order(names(sds))]
    ## TODO: Is this next part even helpful?
    if (all(names(sds) != coefs$rule)) {
      warning("There seems to be a problem with the ordering or size of the
              coefficient and sd vectors. Importances cannot be calculated.")
    }
    
    # baselearner importance is given by abs(coef*st.dev), see F&P section 6):
    if (standardize) {
      baseimps <- data.frame(coefs, sd = sds, imp = abs(coefs$coefficient)*sds/sd_y)
    } else {
      baseimps <- data.frame(coefs, sd = sds, imp = abs(coefs$coefficient)*sds)
    }
    
    
    ## Step 2: Calculate variable importances:
    
    # For factors, importances for each level should be added together.
    # first get indicators for assignments in modmat which are not rules:
    inds <- attr(object$modmat, "assign")[-grep("rule", colnames(object$modmat))]
    # add names in modelframe and modelmatrix to baselearner importances:
    frame.mat.conv <- data.frame(
      modmatname = colnames(object$modmat)[-grep("rule", colnames(object$modmat))],
      modframename = attr(attr(object$data, "terms"), "term.labels")[inds],
      stringsAsFactors = FALSE)
    # We set sort to FALSE to get comparable results across platforms
    baseimps <- base::merge.data.frame(
      frame.mat.conv, baseimps, by.x = "modmatname", by.y = "rule",
      all.x = TRUE, all.y = TRUE, sort = FALSE)
    baseimps <- baseimps[baseimps$coefficient != 0,]
    baseimps <- baseimps[baseimps$description != "(Intercept) ",]
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
    varimps <- data.frame(varname = object$x_names, imp = 0,
                          stringsAsFactors = FALSE)
    # Get importances for rules:
    for(i in 1:nrow(varimps)) { # for every variable:
      # For every baselearner:
      for(j in 1:nrow(baseimps)) {
        # if the variable name appears in the rule:
        #   (Note: EXACT matches are needed, so 1) there should be a space before 
        #     and after the variable name in the rule and thus 2) there should be 
        #     a space added before the description of the rule)
        if(grepl(paste0(" ", varimps$varname[i], " "), paste0(" ", baseimps$description[j]))) {
          # then count the number of times it appears in the rule:
          n_occ <- length(gregexpr(paste0(" ", varimps$varname[i], " "),
                                   paste0(" ", baseimps$description[j]), fixed = TRUE)[[1]])
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
    baseimps <- baseimps[order(baseimps$imp, decreasing = TRUE, method = "radix"),]
    varimps <- varimps[order(varimps$imp, decreasing = TRUE, method = "radix"),]
    varimps <- varimps[varimps$imp != 0,]
    if (plot & nrow(varimps) > 0) {
      if (diag.xlab) {
        xlab.pos <- barplot(height = varimps$imp, xlab = "", ylab = ylab, 
                            main = main, cex.axis = cex.axis, ...)
        ## add specified number of trailing spaces to variable names:
        plotnames <- varimps$varname
        if (diag.xlab.vert > 0) {
          for (i in 1:diag.xlab.vert) {
            plotnames <- paste0(plotnames, " ")
          }
        }
        text(xlab.pos + diag.xlab.hor, par("usr")[3], srt = 45, adj = 1, xpd = TRUE, 
             labels = plotnames, cex = cex.axis)
      } else {
        barplot(height = varimps$imp, names.arg = varimps$varname, ylab = ylab,
                main = main, ...)
      }
    }
    if (!is.na(round)) {
      varimps[,"imp"] <- round(varimps[,"imp"], digits = round)
      baseimps[,c("imp", "coefficient", "sd")] <- round(
        baseimps[,c("imp", "coefficient", "sd")], digits = round)
    }
    row.names(baseimps) <- NULL
    row.names(varimps) <- NULL
    
    return(invisible(list(
      varimps = varimps, 
      baseimps = data.frame(
        rule = baseimps$modmatname,
        baseimps[baseimps$description != "(Intercept) ", c("description", "imp", "coefficient", "sd")],
        stringsAsFactors = FALSE))))
    } else {
      warning("No non-zero terms in the ensemble. All importances are zero.")
      return(invisible(NULL))
    }
}





#' Compute bootstrapped null interaction models
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
#' @seealso \code{\link{pre}}, \code{\link{interact}} 
bsnullinteract <- function(object, nsamp = 10, parallel = FALSE,
                           penalty.par.val = "lambda.1se", verbose = FALSE)
{
  # Preliminaries:
  if(object$family == "binomial") {
    stop("bsnullinteract is not yet available for categorical outcomes.")
  }
  if(parallel) {
    if (!("foreach" %in% installed.packages()[,1])) {
      warning("Parallel computation of function bsnullinteract() requires package foreach,
              which is currently not installed. Argument parallel will be set to FALSE.
              To run in parallel, download and install package foreach from CRAN, and run again.")
      parallel <- FALSE
    }
  }
  # create call for generating bootstrapped null models:
  bsnullmodcall <- object$call
  bsnullmodcall$maxdepth <- 1
  bsnullmodcall$verbose <- verbose
  # create call for model allowing for interactions, grown on bootstrapped
  # datasets without interactions:
  bsintmodcall <- object$call
  bsintmodcall$verbose <- FALSE
  # compute bootstrapped null datasets (i.e., datasets with no interactions):
  if (parallel) {
    if (verbose) cat("This may take a while.")
    bs.ens <- foreach::foreach(i = 1:nsamp) %dopar% {
      # step 1: Take bootstrap sample {x_p, y_p}:
      bs_inds <- sample(1:nrow(object$orig_data), nrow(object$orig_data), replace = TRUE)
      bsdataset <- object$orig_data[bs_inds,]
      # step 2: Build F_A, a null interaction model involving main effects only using {x_p, y_p}:
      bsnullmodcall$data <- bsdataset
      bs.ens.null <- eval(bsnullmodcall)
      # step 3: first part of formula 47 of F&P2008:
      # Calculate predictions F_A(x) for original x, using the null interaction model F_A:
      F_A_of_x <- predict.pre(bs.ens.null, newdata = object$orig_data, 
                              penalty.par.val = penalty.par.val)
      # step 4: third part of formula 47 of F&P2008:
      # Calculate predictions F_A(x_p):
      F_A_of_x_p <- predict.pre(bs.ens.null, newdata = bsdataset,
                                penalty.par.val = penalty.par.val)
      # step 5: Calculate ytilde of formula 47 of F&P2008:
      ytilde <- F_A_of_x + object$data[bs_inds, object$y_name] - F_A_of_x_p
      # step 6: Build a model using (x,ytilde), using the same procedure as was
      # originally applied to (x,y):
      bsintmodcall$data <- object$orig_data
      bsintmodcall$data[,all.vars(object$call$formula[[2]])] <- ytilde
      eval(bsintmodcall)
    }
  } else {
    bs.ens <- list()
    if (verbose) cat("This may take a while. Computing null model ")
    for(i in 1:nsamp) {
      if (verbose) {cat(i, "of", nsamp, ", ")}
      # step 1: Take bootstrap sample {x_p, y_p}:
      bs_inds <- sample(1:nrow(object$orig_data), nrow(object$orig_data), replace = TRUE)
      bsdataset <- object$orig_data[bs_inds,]
      # step 2: Build F_A, a null interaction model involving main effects only using {x_p, y_p}:
      bsnullmodcall$data <- as.symbol(quote(bsdataset))
      bs.ens.null <- eval(bsnullmodcall, envir = environment())
      # step 3: first part of formula 47 of F&P2008:
      # Calculate predictions F_A(x) for original x, using the null interaction model F_A:
      F_A_of_x <- predict.pre(bs.ens.null, newdata = object$orig_data)
      # step 4: third part of formula 47 of F&P2008:
      # Calculate predictions F_A(x_p):
      F_A_of_x_p <- predict.pre(bs.ens.null, newdata = bsdataset,
                                penalty.par.val = penalty.par.val)
      # step 5: Calculate ytilde of formula 47 of F&P2008:
      # FIXME: Does not compute for categorical outcomes: 
      ytilde <- F_A_of_x + object$data[bs_inds, object$y_name] - F_A_of_x_p
      # step 6: Build a model using (x,ytilde), using the same procedure as was
      # originally applied to (x,y):
      tmp <- object$orig_data
      tmp[,all.vars(object$call$formula[[2]])] <- ytilde
      bsintmodcall$data <- as.symbol(quote(tmp))
      bs.ens[[i]] <- eval(bsintmodcall, envir = environment())
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
#' @param nullmods object with bootstrapped null interaction models, resulting
#' from application of \code{bsnullinteract}.
#' @param penalty.par.val character. Which value of the penalty parameter
#' criterion should be used? The value yielding minimum cv error
#' (\code{"lambda.min"}) or penalty parameter yielding error within 1 standard
#' error of minimum cv error ("\code{lambda.1se}")? Alternatively, a numeric 
#' value may be specified, corresponding to one of the values of lambda in the 
#' sequence used by glmnet, for which estimated cv error can be inspected by 
#' running \code{object$glmnet.fit} and \code{plot(object$glmnet.fit)}.
#' @param quantprobs numeric vector of length two. Probabilities that should be
#' used for plotting the range of bootstrapped null interaction model statistics.
#' Only used when \code{nullmods} argument is specified and \code{plot = TRUE}.
#' The default yields sample quantiles corresponding to .05 and .95 probabilities.  
#' @param plot logical. Should interaction statistics be plotted?
#' @param col character vector of length one or two. Color for plotting 
#' interaction statistics. The first color specified is used to plot the 
#' interaction statistic from the training data, the second color specifed
#' is used to plot the interaction statistic distribution from the bootstrapped
#' null interaction models. Only used when \code{plot = TRUE}. Only the first 
#' element of vector is used if \code{nullmods = NULL}.
#' @param ylab character string. Label to be used for plotting y-axis.
#' @param main character. Main title for the bar plot.
#' @param  se.linewidth numeric. Width of the whiskers of the plotted standard 
#' error bars (in inches).
#' @param k integer. Calculating interaction test statistics is a computationally
#' intensive, so  calculations are split up in several parts to prevent memory
#' allocation errors. If a memory allocation error still occurs, increase k.
#' @param verbose logical. Should progress information be printed to the
#' command line?
#' @param parallel logical. Should parallel foreach be used? Must register
#' parallel beforehand, such as doMC or others.
#' @param ... Additional arguments to be passed to \code{barplot}.
#' @examples
#' \donttest{
#'  set.seed(42)
#'  airq.ens <- pre(Ozone ~ ., data=airquality[complete.cases(airquality),])
#'  interact(airq.ens, c("Temp", "Wind", "Solar.R"))}
#' @details Can be computationally intensive, especially when nullmods is 
#' specified, in which case setting \verb{parallel = TRUE} may improve speed.
#' @return Function \code{interact()} returns and plots interaction statistics
#' for the specified predictor variables. If nullmods is not specified, it 
#' returns and plots only the interaction test statistics for the specified 
#' fitted prediction rule ensemble. If nullmods is specified, the function 
#' returns a list, with elements \code{$fittedH2}, containing the interaction
#' statistics of the fitted ensemble, and \code{$nullH2}, which contains the
#' interaction test statistics for each of the bootstrapped null interaction 
#' models.
#'  
#' If \code{plot = TRUE} (the default), a barplot is created with the 
#' interaction test statistic from the fitted prediction rule ensemble. If 
#' \code{nullmods} is specified, bars representing the median of the 
#' distribution of interaction test statistics of the bootstrapped null 
#' interaction models are plotted. In addition, error bars representing the
#' quantiles of the distribution (their value specified by the \code{quantprobs} 
#' argument) are plotted. These allow for testing the null hypothesis of no 
#' interaction effect for each of the input variables. 
#' 
#' Note that the error rates of null hypothesis tests of interaction effects 
#' have not yet been studied in detail, but likely depend on the number of 
#' generated bootstrapped null interaction models as well as the complexity of 
#' the fitted ensembles. Users are therefore advised to test for the presence 
#' of interaction effects by setting the \code{nsamp} argument of the function 
#' \code{bsnullinteract} \eqn{\geq 100} (even though this may take a lot of 
#' computation time). Also, users are advised to test for the presence of 
#' interactions only with fitted ensembles that are neither too sparse nor too 
#' complex, that is, ensembles that are selected by setting the 
#' \code{penalty.par.val} argument equal to \code{"lambda.min"} or 
#' \code{"lambda.1se"}.  
#' @export
#' @seealso \code{\link{pre}}, \code{\link{bsnullinteract}} 
interact <- function(object, varnames = NULL, nullmods = NULL, 
                     penalty.par.val = "lambda.1se", quantprobs = c(.05, .95),
                     plot = TRUE, col = c("yellow", "blue"), 
                     ylab = "Interaction strength", 
                     main = "Interaction test statistics", 
                     se.linewidth = .05,
                     parallel = FALSE, k = 10, verbose = FALSE, ...) {
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
      medians <- rbind(H, apply(nullH, 2, mean))
      H0_medians <- apply(nullH, 2, median)
      lower_quant <- apply(nullH, 2, quantile, probs = quantprobs[1])
      upper_quant <- apply(nullH, 2, quantile, probs = quantprobs[2])
      x_coords <- barplot(medians, beside = TRUE, 
                          ylim = c(0, max(upper_quant, medians)), 
                          las = 1, main = main, col = col, ...)
      x_coords <- x_coords[!1:nrow(x_coords)%%2,] 
      segments(x_coords, lower_quant, x_coords, upper_quant)
      arrows(x_coords, lower_quant, x_coords, upper_quant, lwd = 1.5, angle = 90, 
             code = 3, length = se.linewidth)
    }
  }
  if(is.null(nullmods)) {
    return(H)
  } else {
    return(list(fittedH2 = H, nullH2 = nullH))
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
#' \code{linear.terms = FALSE}) being plotted. Default is \code{NULL}, 
#' resulting in all terms of the final ensemble to be plotted.
#' @param plot.dim integer vector of length two. Specifies the number of rows
#' and columns in the plot. The default yields a plot with three rows and three 
#' columns, depicting nine baselearners per plot. If 
#' \code{nterms > plot.dim[1] * plot.dim[2]}, multiple plotting pages will be 
#' created.
#' @param ask logical. Should user be prompted before starting a new page of
#' plots?
#' @param exit.label character string. Label to be printed in nodes to which 
#' the rule does not apply (``exit nodes'')?
#' @param standardize logical. Should printed importances be standardized? See
#' \code{\link{importance}}.
#' @param ... Arguments to be passed to \code{\link[grid]{gpar}}.
#' @examples
#' \donttest{
#'  set.seed(42)
#'  airq.ens <- pre(Ozone ~ ., data = airquality[complete.cases(airquality),])
#'  plot(airq.ens)}
#' @export
#' @seealso \code{\link{pre}}, \code{\link{print.pre}}
#' @method plot pre
plot.pre <- function(x, penalty.par.val = "lambda.1se", linear.terms = TRUE, 
                     nterms = NULL, ask = FALSE, exit.label = "0", 
                     standardize = FALSE, plot.dim = c(3, 3), ...) {
  ## Preliminaries:
  if (!("grid" %in% installed.packages()[,1])) {
    stop("Function plot.pre requires package grid. Download and install package
         grid from CRAN, and run again.")
  }
  
  
  
  x <- airq.ens
  penalty.par.val = "lambda.1se"
  linear.terms = FALSE 
  nterms = 6
  ask = FALSE
  exit.label = "0" 
  standardize = FALSE
  plot.dim = c(3, 3)
  
  
  ## Get nonzero terms:
  nonzeroterms <- importance(x, plot = FALSE, global = TRUE, 
                             penalty.par.val = penalty.par.val, 
                             standardize = standardize)$baseimps
  if (!linear.terms) {
    nonzeroterms <- nonzeroterms[grep("rule", nonzeroterms$rule),]
  }
  if (!is.null(nterms)) {
    nonzeroterms <- nonzeroterms[1:nterms,]
  }

  ## Grab baselearner components:
  conditions <- list()
  for(i in 1:nrow(nonzeroterms)) { # i is a counter for terms
    if (length(grep("&", nonzeroterms$description[i], )) > 0) { # get rules with multiple conditions:
      conditions[[i]] <- unlist(strsplit(nonzeroterms$description[i], split = " & "))
    } else if (!grepl("rule", nonzeroterms$rule[i])) {
      conditions[[i]] <- "linear" # flags linear terms
    } else {
      conditions[[i]] <- nonzeroterms$description[i] # gets rules with only one condition
    }
  }
  
  ## for every non-zero term, calculate the number of the plot, row and column where it should appear.
  n_terms_per_plot <- plot.dim[1] * plot.dim[2]
  nplots <- ceiling(nrow(nonzeroterms) / n_terms_per_plot)
  nonzeroterms$plotno <- rep(1:nplots, each = n_terms_per_plot)[1:nrow(nonzeroterms)]
  nonzeroterms$rowno <- rep(rep(1:plot.dim[1], each = plot.dim[2]), length.out = nrow(nonzeroterms))
  nonzeroterms$colno <- rep(rep(1:plot.dim[2], times = plot.dim[1]), length.out = nrow(nonzeroterms))
  
  ## Generate a plot for every baselearner:
  for(i in 1:nrow(nonzeroterms)) {
    if (conditions[[i]][1] == "linear") { 
      
      ## Plot linear term:
      
      ## Open new plotting page if needed:
      if (nonzeroterms$rowno[i] == 1 && nonzeroterms$rowno[i] == 1) {
        grid::grid.newpage()
        grid::pushViewport(grid::viewport(layout = grid::grid.layout(plot.dim[1], plot.dim[2])))
      }
      ## open correct viewport:
      grid::pushViewport(grid::viewport(layout.pos.col = nonzeroterms$colno[i],
                                        layout.pos.row = nonzeroterms$rowno[i]))
      ## Plot the linear term:
      grid::grid.text(paste0("Linear effect of ", nonzeroterms$rule[i], 
                             "\n\n Coefficient = ", round(nonzeroterms$coefficient[i], digits = 3),
                             "\n\n Importance = ", round(nonzeroterms$imp[i], digits = 3)),
                      gp = grid::gpar(...))
      grid::popViewport()
      
    } else { 
      
      ## Otherwise, plot rule:
      
      # Create lists of arguments and operators for every condition:
      cond <- list()
      # check whether the operator is " < ", " <= " or "%in%  "
      # split the string using the operator, into the variable name and splitting value, 
      # which is used to define split = partysplit(id, value)
      # make it a list:
      for (j in 1:length(conditions[[i]])) {
        condition_j <- conditions[[i]][[j]]
        cond[[j]] <- character()
        if (length(grep(" > ", condition_j)) > 0) {
          cond[[j]][1] <- unlist(strsplit(condition_j, " > "))[1]
          cond[[j]][2] <- " > "
          cond[[j]][3] <- unlist(strsplit(condition_j, " > "))[2]
        }
        if (length(grep(" <= ", condition_j)) > 0) {
          cond[[j]][1] <- unlist(strsplit(condition_j, " <= "))[1]
          cond[[j]][2] <- " <= "
          cond[[j]][3] <- unlist(strsplit(condition_j, " <= "))[2]
        }
        if (length(grep(" %in% ", condition_j)) > 0) {
          cond[[j]][1] <- unlist(strsplit(condition_j, " %in% "))[1]
          cond[[j]][2] <- " %in% "
          cond[[j]][3] <- unlist(strsplit(condition_j, " %in% "))[2]
        }
      }
      ncond <- length(cond)
      cond <- rev(cond)
      
      ## generate empty dataset for all the variables appearing in the rules: 
      treeplotdata <- data.frame(matrix(ncol = ncond))
      for (j in 1:ncond) {
        names(treeplotdata)[j] <- cond[[j]][1]
        if (cond[[j]][2] == " %in% ") {
          treeplotdata[,j] <- factor(treeplotdata[,j])
          faclevels <- substring(cond[[j]][3], first = 2)
          faclevels <- gsub(pattern = "\"", replacement = "", x = faclevels, fixed = TRUE)
          faclevels <- gsub(pattern = "(", replacement = "", x = faclevels, fixed = TRUE)
          faclevels <- gsub(pattern = ")", replacement = "", x = faclevels, fixed = TRUE)
          faclevels <- unlist(strsplit(faclevels, ", ",))
          levels(treeplotdata[,j]) <- c(
            levels(x$data[,cond[[j]][1]])[levels(x$data[,cond[[j]][1]]) %in% faclevels],
            levels(x$data[,cond[[j]][1]])[!(levels(x$data[,cond[[j]][1]]) %in% faclevels)])
          cond[[j]][3] <- length(faclevels)
        }
      }

      
      ## Generate partynode objects for plotting:
      
      nodes <- list()
      
      ## Create level 0 (bottom level, last two nodes):
      ## If last condition has " > " : exit node on left, coefficient on right:
      if (cond[[1]][2] == " > ") { # If condition involves " > ", the tree has nonzero coef on right:
        nodes[[1]] <- list(id = 1L, split = NULL, kids = NULL, surrogates = NULL, 
                           info = exit.label)
        nodes[[2]] <- list(id = 2L, split = NULL, kids = NULL, surrogates = NULL,
                           info = round(nonzeroterms$coefficient[i], digits = 3))
      } else { 
      ## If last condition has " <= " or " %in% " : coefficient on left, exit node on right:
        nodes[[1]] <- list(id = 1L, split = NULL, kids = NULL, surrogates = NULL,
                           info = round(nonzeroterms$coefficient[i], digits = 3))
        nodes[[2]] <- list(id = 2L, split = NULL, kids = NULL, surrogates = NULL,
                           info = exit.label)
      }
      class(nodes[[1]]) <- class(nodes[[2]]) <- "partynode"
      
      
      ## Create inner levels (if necessary):
      if (ncond > 1) {
        for (level in 1:(ncond - 1)) {
          if (cond[[level + 1]][2] == " > ") { 
            ## If condition in level above has " > " : exit node on left, right node has kids:
            nodes[[level * 2 + 1]] <- list(id = as.integer(level * 2 + 1), 
                                          split = NULL, 
                                          kids = NULL, 
                                          surrogates = NULL, 
                                          info = exit.label)
            nodes[[level * 2 + 2]] <- list(id = as.integer(level * 2 + 2), 
                                          split = partysplit(as.integer(level), breaks = as.numeric(cond[[level]][3])),
                                          kids = list(nodes[[level * 2 - 1]], nodes[[level * 2]]),
                                          surrogates = NULL, 
                                          info = NULL)
            } else { 
            ## If condition in level above has " <= " or " %in% " : left node has kids, exit node right:
            nodes[[level * 2 + 1]] <- list(id = as.integer(level * 2 + 1),
                                          split = partysplit(as.integer(level), breaks = as.numeric(cond[[level]][3])),
                                          kids = list(nodes[[level * 2 - 1]], nodes[[level * 2]]),
                                          surrogates = NULL, 
                                          info = NULL)
            nodes[[level * 2 + 2]] <- list(id = as.integer(level * 2 + 2), 
                                          split = NULL,
                                          kids = NULL, 
                                          surrogates = NULL, 
                                          info = exit.label)
         }  
         class(nodes[[level * 2 + 1]]) <- class(nodes[[level * 2 + 2]]) <- "partynode"
      }
    }
      
    ## Create root node:
    nodes[[ncond * 2 + 1]] <- list(id = as.integer(ncond * 2 + 1),
                                   split = partysplit(as.integer(ncond), breaks = as.numeric(cond[[ncond]][3])),
                                   kids = list(nodes[[ncond * 2 - 1]], nodes[[ncond * 2]]),
                                   surrogates = NULL, 
                                   info = NULL)
    class(nodes[[ncond * 2 + 1]]) <- "partynode"
      
    ## Open new plotting page if needed:
    if (nonzeroterms$rowno[i] == 1 && nonzeroterms$colno[i] == 1) {
      grid::grid.newpage()
      grid::pushViewport(grid::viewport(layout = grid::grid.layout(plot.dim[1], plot.dim[2])))
    }
    
    ## Open viewport:
    grid::pushViewport(grid::viewport(layout.pos.col = nonzeroterms$colno[i],
                                      layout.pos.row = nonzeroterms$rowno[i]))
    
    ## pPlot the rule:
    fftree <- party(nodes[[ncond * 2 + 1]], data = treeplotdata)
    plot(fftree, newpage = FALSE, 
         main = paste0(nonzeroterms$rule[i], ": Importance = ", round(nonzeroterms$imp[i], digits = 3)),
         inner_panel = node_inner(fftree, id = FALSE),
         terminal_panel = node_terminal(fftree, id = FALSE))#, gp = grid::gpar(...))
    grid::popViewport()
  }
  }
  
  if (ask) {
    grDevices::devAskNewPage(ask = FALSE)
  }
}



## Internal function for creating legend for corplot:
image.scale <- function(z, col, breaks, axis.pos = 4, add.axis = TRUE) {
  poly <- vector(mode = "list", length(col))
  for(i in seq(poly)) {
    poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
  }
  if (axis.pos %in% c(1,3)) { 
    ylim <- c(0,1)
    xlim <- range(breaks)
  }
  if (axis.pos %in% c(2,4)){
    ylim <- range(breaks)
    xlim <- c(0,1)
  }
  plot(1, 1, t = "n", ylim = ylim, xlim = xlim, axes = FALSE, xlab = "", 
       ylab = "", xaxs = "i", yaxs = "i")  
  for(i in seq(poly)){
    if(axis.pos %in% c(1,3)){
      polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
    }
    if(axis.pos %in% c(2,4)){
      polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
    }
  }
}



#' Plotting baselearner correlations
#' 
#' \code{corplot} plots correlations between baselearners
#'  
#' @param object object of class pre
#' @param penalty.par.val character. Value of the penalty parameter value 
#' \eqn{\lambda} to be used for selecting the final ensemble. The ensemble 
#' with penalty parameter criterion yielding minimum cv error 
#' (\code{"lambda.min"}) is taken, by default. Alternatively, the penalty 
#' parameter yielding error within 1 standard error of minimum cv error 
#' ("\code{lambda.1se}"), or a numeric value may be specified, corresponding 
#' to one of the values of lambda in the sequence used by glmnet,
#' for which estimated cv error can be inspected by running \code{x$glmnet.fit}
#' and \code{plot(x$glmnet.fit)}.
#' @param colors vector of contiguous colors to be used for plotting. If 
#' \code{colors = NULL} (default), \code{colorRampPalette(c("#053061", "#2166AC", 
#' "#4393C3", "#92C5DE", "#D1E5F0", "#FFFFFF", "#FDDBC7", "#F4A582", "#D6604D", 
#' "#B2182B", "#67001F"))(200)} is used. A different set of plotting colors can 
#' be specified, for example: \code{colors = cm.colors(100)}, or
#' \code{colorRampPalette(c("blue", "white", "red"))(150)}. See
#' \code{\link[grDevices]{cm.colors}} or \code{\link[grDevices]{colorRampPalette}}.
#' @param fig.plot plotting region to be used for correlation plot. See 
#' \code{fig} under \code{\link{par}}.
#' @param fig.legend plotting region to be used for legend. See \code{fig} 
#' under \code{\link{par}}.
#' @param legend.breaks numeric vector of breakspoints and colors to be 
#' depicted in the plot's legend.
#' @examples \donttest{
#' set.seed(42)
#' airq.ens <- pre(Ozone ~ ., data = airquality[complete.cases(airquality),])
#' corplot(airq.ens)
#' }
#' @export
corplot <- function(object, penalty.par.val = "lambda.1se", colors = NULL,
                    fig.plot = c(0, 0.85, 0, 1), fig.legend = c(.8, .95, 0, 1),
                    legend.breaks = seq(-1, 1, by = .1)) {
  
  if (is.null(colors)) {
    colors <- grDevices::colorRampPalette(
      c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#FFFFFF", 
        "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"))(200)
  }
  ## get coefficients:
  coefs <- as(coef.glmnet(object$glmnet.fit, s = penalty.par.val), 
              Class = "matrix")
  ## create correlation matrix of non-zero coefficient learners:
  cormat <- cor(object$modmat[,names(coefs[coefs != 0,])[-1]])
  ## set layout for plotting correlation matrix and legend::
  layout(matrix(rep(c(1, 1, 1, 1, 2), times = 5), 5, 5, byrow = TRUE))
  ## set region for correlation matrix:
  par(fig = fig.plot)
  ## plot correlation matrix:
  image(x = 1:nrow(cormat), y = 1:ncol(cormat), z = unlist(cormat), 
        axes = FALSE, zlim = c(-1, 1),  
        xlab = "", ylab = "", srt = 45, 
        #col = grDevices::cm.colors(21),
        col = colors)
  axis(1, at = 1:nrow(cormat), labels = colnames(cormat), las = 2)
  axis(2, at = 1:ncol(cormat), labels = colnames(cormat), las = 2)
  ## plot legend:
  par(fig = fig.legend, new = TRUE)
  col_inds <- round(seq(1, length(colors), length.out = length(legend.breaks)-1))
  image.scale(z = legend.breaks, 
              #col = grDevices::cm.colors(11), 
              col = colors[col_inds],
              axis.pos = 4,
              breaks = legend.breaks)
  axis(4, at = legend.breaks, las = 2)
  return(invisible(cormat))
}

