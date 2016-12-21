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
#' @param sampfrac numeric value greater than 0, and smaller than or equal to 1. 
#' Fraction of randomly selected training observations used to produce each tree. 
#' Setting this to values < 1 will result in subsamples being drawn without 
#' replacement (i.e., subsampling). Setting this equal to 1 will result in 
#' bootstrap sampling.
#' @param seed numeric. Random seed to be used in deriving the final ensemble 
#' (for reproducability).
#' @param maxdepth numeric. Maximal depth of trees to be grown. Defaults to 3,
#' resulting in trees with max 15 nodes (8 terminal and 7 inner nodes), and 
#' therefore max 15 rules.
#' @param learnrate numeric. Learning rate for sequentially induced trees.
#' @param removeduplicates logical. Remove rules from the ensemble which have 
#' the exact same support in training data?
#' @param maxrules numeric. Approximate maximum number of rules to be generated. 
#' The number of rules in the final ensemble will be smaller, due to the omission 
#' of rules with identical conditions or support.
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
#' deriving the final ensemble. The default is \code{type.measure = "deviance"}, 
#' which uses squared-error for gaussian models (a.k.a. \code{type.measure = 
#' "mse"}). \code{type.measure = "mse"} or \code{type.measure = "mae"} (mean 
#' absolute error) measure the deviation from the fitted mean to the response.
#' @param verbose logical. Should information on the initial and final ensemble 
#' be printed to the command line?
#' @param ctreecontrol A list with control parameters, see 
#' \code{link[partykit]{ctree_control}}. Ignored for nominal output variables 
#' when \code{learnrate} > 0.
#' @param ... Additional arguments to be passed to 
#' \code{\link[glmnet]{cv.glmnet}}.
#' @return an object of class \code{pre}, which is a list with many elements 
#' @details Inputs can be continuous, ordered or factor variables. Continuous 
#' variables
#' @examples \donttest{
#' airq.ens <- pre(Ozone ~ ., data = airquality[complete.cases(airquality),])}
#' @import glmnet partykit datasets
#' @export
pre <- function(formula, data, type = "both", weights = rep(1, times = nrow(data)), 
                 sampfrac = .5, seed = 42, maxdepth = 3, learnrate = 0.01, 
                 removeduplicates = TRUE, maxrules = 2000, mtry = Inf, 
                 thres = 1e-07, standardize = FALSE, winsfrac = .025, 
                 normalize = TRUE, nfolds = 10, mod.sel.crit = "deviance", 
                 verbose = TRUE, ctreecontrol = ctree_control(),
                 ...)   
{
  ###################
  ## Preliminaries ##
  ###################
  
  if(!is.data.frame(data)) {
    stop("data should be a data frame.")
  }
  if(length(sampfrac) != 1 || sampfrac < 0.01 || sampfrac > 1) {
    stop("Bad value for 'sampfrac'")
  }
  if(length(type) != 1 || (type != "rules" & type != "both" & type != "linear")) {
    stop("Argument type should equal 'both', 'rules' or 'linear'")
  }
  if(length(winsfrac) != 1 || winsfrac < 0 || winsfrac > 0.5) {
    stop("Bad value for 'winsfrac'.")
  }
  if(!is.logical(verbose)) {
    stop("Bad value for 'verbose'.")
  }  
  set.seed(seed)
  data <- model.frame(formula, data, na.action = NULL)
  x_names <- attr(attr(data, "terms"), "term.labels")
  y_name <- names(data)[!names(data)%in%x_names]
  if(is.factor(data[,y_name])) {
    classify <- TRUE
  } else {
    classify <- FALSE
  }
  n <- nrow(data)
  if(!(is.numeric(data[,y_name]) | is.factor(data[,y_name]))) {
    stop("Output variable should be continuous (class numeric) or binary (class 
         factor)")
  }
  if(nlevels(data[,y_name]) > 2) {
    stop("No support for multinomial output variables yet.")
  }
  if(any(sapply(data[,x_names], is.character))) {
    stop("Variables specified in formula and data argument are of class character. 
         Please coerce to class 'numeric', 'factor' or 'ordered' 'factor':", 
         x_names[sapply(data[,x_names], is.character)])
  }
  if(classify & learnrate != 0 & !is.infinite(mtry)) {
    warning("Value specified for mtry will not be used when the outcome variable
            is binary and learnrate > 0")
  }
  if(any(is.na(data))) {
    data <- data[complete.cases(data),]
    n <- nrow(data)
    warning("Some observations have missing values and will be removed. 
            New sample size is ", n, ".\n", immediate. = TRUE)
  }
  if(verbose) {
    if(classify) {
      cat("A rule ensemble for prediction of a binary output variable will be 
          created.\n")
    } else {
      cat("A rule ensemble for prediction of a continuous output variable will 
          be created.\n")
    }
  }  
  
  #############################
  ## Derive prediction rules ##
  #############################
  
  if(type != "linear") {
    rules <- vector()
    treecount <- 0
    # If learning rate = 0, always use ctree():
    if(learnrate == 0) {
      while(length(rules) <= maxrules) {
        # Take subsample of dataset
        treecount <- treecount + 1
        if(sampfrac == 1) { # then bootstrap
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
        rules <- append(rules, unlist(partykit:::.list.rules.party(tree)))
      }
    }
    if(learnrate > 0) {
      if(!classify) {
        y_learn <- data[,y_name]
        while(length(rules) <= maxrules) {
          # Take subsample of dataset
          treecount <- treecount + 1
          if(sampfrac == 1) { # then bootstrap
            subsample <- sample(1:n, size = n, replace = TRUE)
          } else { # else subsample
            subsample <- sample(1:n, size = round(sampfrac * n), replace = FALSE)
          }
          subsampledata <- data[subsample,]
          subsampledata[,y_name] <- y_learn[subsample]
          # Make sure ctree() can find object specified by weights argument: 
          environment(formula) <- environment()
          # Grow ctree on subsample:
          tree <- ctree(formula, data = subsampledata, weights = weights[subsample], 
                        maxdepth = maxdepth, mtry = mtry)
          # Collect rules from tree:
          rules <- append(rules, unlist(partykit:::.list.rules.party(tree)))
          # Substract predictions from current y:
          y_learn <- y_learn - learnrate * predict(tree, newdata = data)
        }
      }
      if(classify) { # if ouput var is a factor:
        data2 <- data.frame(data, offset = 0)
        while(length(rules) <= maxrules) {
          # Take subsample of dataset:
          treecount <- treecount + 1
          if(sampfrac == 1) { # then bootstrap:
            subsample <- sample(1:n, size = n, replace = TRUE)
          } else { # else subsample:
            subsample <- sample(1:n, size = round(sampfrac * n), replace = FALSE)
          }
          subsampledata <- data2[subsample,]
          # Make sure ctree() can find object specified by weights argument: 
          environment(formula) <- environment()
          # Grow ctree on subsample:
          tree <- glmtree(formula(paste(paste(y_name, " ~ 1 |"), 
                                        paste(x_names, collapse = "+"))),
                          data = subsampledata, family = "binomial",
                          weights = weights[subsample], maxdepth = maxdepth+1,  
                          offset = offset)
          # Collect rules from tree:
          rules <- append(rules, unlist(partykit:::.list.rules.party(tree)))
          # Update offset:
          data2$offset <- data2$offset + learnrate * predict(
            tree, newdata = data2, type = "link")
        }
      } 
    }
    nrules <- length(rules)
    if(verbose){
      cat("\nA total of", treecount, "trees were grown, and a total of", nrules, 
          "rules were generated initially.")
    }
    # Keep unique, non-empty rules only:
    rules <- unique(rules[!rules==""])
    if(verbose) {
      cat("\n\nA total of", nrules - length(rules), "rules were duplicate or empty
            rules and removed from the initial ensemble.")
    }
    # Create dataframe with 0-1 coded rules:
    if(length(rules)>0) {
      rulevars <- data.frame(
        rule1 = as.numeric(with(data, eval(parse(text = rules[[1]])))))
      for(i in 2:length(rules)) {
        rulevars[,paste("rule", i, sep="")] <- as.numeric(
          with(data, eval(parse(text = rules[[i]]))))
      }
      if(removeduplicates) {
        # Remove rules with identical support:
        duplicates <- duplicated(t(rulevars))
        duplicates.removed <- data.frame(name = colnames(rulevars)[duplicates],
                                         description = rules[duplicates])
        rulevars <- rulevars[,!duplicates]
        rules <- rules[!duplicates]
        if(verbose) {
          cat("\n\nA total of", sum(duplicates), "generated rules had 
              support identical to earlier rules and were removed from the initial 
              ensemble ($duplicates.removed shows which, if any).")
        }
      } else {
        duplicates.removed <- NULL
      }
      if(verbose) {
        cat("\n\nAn initial ensemble consisting of", ncol(rulevars), "rules was 
            succesfully created.")  
      }
    }
  }
  if(length(rules)==0) {
    type <- "linear"
    warning("No prediction rules could be derived from dataset. Final ensemble 
            will include linear terms only, type will be set to 'linear'.")
  }
  
  ######################################################
  ## Prepare rules, linear terms and outcome variable ##
  ######################################################
  
  x <- data[,x_names]

  # convert ordered categorical predictor variables to linear terms:
  x[,sapply(x, is.ordered)] <- as.numeric(as.character(x[,sapply(x, is.ordered)]))

  if(type == "rules") {
    x <- rulevars
    x_scales <- NULL
  } else { # if type is not rules, linear terms should be prepared:
    # Winsorize numeric variables (section 5 of F&P(2008)):
    if(winsfrac > 0) {
      for(i in 1:ncol(x)) {
        if(is.numeric(x[,i])) { 
          lim <- quantile(x[,i], probs = c(winsfrac, 1-winsfrac))
          x[x[,i] < lim[1],i] <- lim[1]
          x[x[,i] > lim[2],i] <- lim[2]
        }
      }
    }
    # normalize numeric variables:
    if(normalize) { 
      # Normalize linear terms (section 5 of F&P08), if there are any:
      if(sum(sapply(x, is.numeric)) > 0) {
        x_scales <- sapply(x[,sapply(x, is.numeric)], sd, na.rm = TRUE)/0.4
        x[,sapply(x, is.numeric)] <- scale(x[,sapply(x, is.numeric)], 
                                           center = FALSE, scale = x_scales)
      } else {
        x_scales <- NULL
      }
    } else {
      x_scales <- NULL
    } 
    # If both rules and linear terms are in ensemble, combine both:
    if(type == "both") {
      x <- data.frame(x, rulevars)
    }
  }
  modmat.formula <- formula(
    paste(paste(y_name, " ~ -1 +"), paste(colnames(x), collapse = "+")))
  x <- model.matrix(modmat.formula, data = data.frame(data[y_name], x))
  y <- data[,y_name]
  
  ##################################################
  ## Perform penalized regression on the ensemble ##
  ##################################################
  
  if(classify) {
    family <- "binomial"
  } else {
    family <- "gaussian"
  }
  
  glmnet.fit <- cv.glmnet(x, y, nfolds = nfolds, standardize = standardize, 
                          type.measure = mod.sel.crit, thres = thres, 
                          weights = weights, family = family, ...)
  
  ####################
  ## Return results ##
  ####################
  
  lmin_ind <- which(glmnet.fit$lambda == glmnet.fit$lambda.min)
  l1se_ind <- which(glmnet.fit$lambda == glmnet.fit$lambda.1se)
  if(verbose) {
    cat("\n\nFinal ensemble with minimum cv error: \n  lambda = ", 
        glmnet.fit$lambda[lmin_ind], "\n  number of terms = ", 
        glmnet.fit$nzero[lmin_ind], "\n  mean cv error (se) = ", 
        glmnet.fit$cvm[lmin_ind], " (", glmnet.fit$cvsd[lmin_ind], ")", 
        "\n\nEnsemble with cv error within 1se of minimum: \n  lambda = ", 
        glmnet.fit$lambda[l1se_ind],  "\n  number of terms = ", 
        glmnet.fit$nzero[l1se_ind], "\n  mean cv error (se) = ", 
        glmnet.fit$cvm[l1se_ind], " (", glmnet.fit$cvsd[l1se_ind], ")\n", sep="")
  }
  result <- list(glmnet.fit = glmnet.fit, call = match.call(), weights = weights, 
                 data = data, normalize = normalize, x_scales = x_scales, 
                 type = type, x_names = x_names, y_name = y_name, 
                 winsfrac = winsfrac, modmat = x, classify = classify, 
                 modmat.formula = modmat.formula)
  if(type != "linear") {
    result$duplicates.removed <- duplicates.removed
    result$rules <- data.frame(rule = names(rulevars), description = rules)
    result$rulevars <- rulevars 
  }
  class(result) <- "pre"
  return(result)
}



#' Full k-fold cross validation of a pre
#' 
#' \code{cvpre} performs k-fold cross validation on the dataset used to create 
#' the ensemble, providing an estimate of predictive accuracy on future observations.
#' 
#' @param object An object of class \code{\link{pre}}.
#' @param k integer. The number of cross validation folds to be used.
#' @param seed integer. Random seed to be used for generating 
#' cross-validation subsamples.
#' @param verbose logical. Should progress of the cross validation be printed 
#' to the command line?
#' @param pclass numeric. Only used for classification. Cut-off value between 
#' 0 and 1 to be used for classifying to second class. 
#' @return A list with three elements: cvpreds (a vector with cross-validated
#' predicted y values), ss (a vector indicating the cross-validation subsample 
#' each training observation was assigned to) and accuracy. For continuous 
#' outputs, accuracy is a list with elements MSE (mean squared error on test 
#' observations), MAE (mean absolute error on test observations) and 
#' cor_true_pred_y, which is the correlation between the actual and predicted 
#' outcome variable values. For classification, accuracy is a list with elements 
#' SEL (mean squared error on predicted probabilities), AEL (mean absolute 
#' error on predicted probabilities), MCR (average misclassification error rate) 
#' and table (a table with proportions of (in)correctly classified observations 
#' per class).
#' @examples \donttest{
#' airq.ens <- pre(Ozone ~ ., data = airquality[complete.cases(airquality),])
#' airq.cv <- cvpre(airq.ens)}
cvpre <- function(object, k = 10, seed = 42, verbose = TRUE, pclass = .5) {
  set.seed(seed)
  ss <- sample(rep(1:k, length.out = nrow(object$data)), size = nrow(object$data), replace = FALSE)
  #ss <- sample(1:k, size = nrow(object$data), replace = TRUE)   
  cvobject <- list()
  cvpreds <- rep(NA, times = nrow(object$data))
  if(verbose) {
    cat("Running cross validation in fold ")
  }
  for (i in 1:k){
    if(verbose) {
      cat(i, "of", k, ",")
    }
    cl <- object$call
    cl$verbose <- FALSE
    cl$data <- object$data[ss!=i,]
    cvobject[[i]] <- eval(cl)
    if(object$classify) {
      cvpreds[ss==i] <- predict.pre(cvobject[[i]], newdata = object$data[ss==i,], 
                                     type = "response")
    } else {
        cvpreds[ss==i] <- predict.pre(cvobject[[i]], newdata = object$data[ss==i,])
    }
    if(verbose & i==k) {
      cat(" done!\n")
    }
  }
  accuracy <- list()
  if(object$classify) {
    accuracy$SEL<- mean((as.numeric(object$data[,object$y_name]) - 1 - cvpreds)^2)
    accuracy$AEL <- mean(abs(as.numeric(object$data[,object$y_name]) - 1 - cvpreds))
    cvpredsd <- as.numeric(cvpreds > .5)
    accuracy$MCR <- 1 - sum(diag(prop.table(table(cvpredsd, 
                                                  object$data[,object$y_name]))))
    accuracy$table <- prop.table(table(cvpredsd, object$data[,object$y_name]))
  }
  else {
    accuracy$MSE <- mean((object$data[,object$y_name] - cvpreds)^2)
    accuracy$MAE <- mean(abs(object$data[,object$y_name] - cvpreds))
    accuracy$cor_true_pred_y <- cor(cvpreds, object$data[,object$y_name])
  }
  result <- list(cvpreds = cvpreds, ss_indicators = ss, accuracy = accuracy)
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
#' ("\code{lambda.1se}").
#' @param print logical. Should coefficients of the base learners with non-zero 
#' coeffcients in the final ensemble be printed to the command line?
#' @param ... additional arguments to be passed to \code{\link[glmnet]{coef.glmnet}}.
#' @return returns a dataframe with 3 columns: coefficients, rule (rule or 
#' variable name) and description (\code{NA} for linear terms, conditions for 
#' rules). In the command line, the non zero coefficients are printed (when 
#' \code{print = TRUE}).
#' @examples \donttest{
#' airq.ens <- pre(Ozone ~ ., data=airquality[complete.cases(airquality),])
#' coefs <- coef(airq.ens)}
#' @export
#' @method coef pre
coef.pre <- function(object, penalty.par.val = "lambda.1se", print = TRUE, ...)
{
  coefs <- as(coef.glmnet(object$glmnet.fit, s = penalty.par.val, ...), 
              Class = "matrix")
  # coefficients for normalized variables should be unnormalized: 
  if(object$normalize & !is.null(object$x_scales) & object$type != "rules") {
    coefs[names(object$x_scales),] <- coefs[names(object$x_scales),] /
      object$x_scales
  }
  coefs <- data.frame(coefficient = coefs[,1], rule = rownames(coefs))
  if(object$type != "linear") {
    coefs <- merge(coefs, object$rules, all.x=T)
  }
  if(print) {
    nonzerocoefs <- coefs[coefs$coefficient != 0,]
    nonzerocoefs$coefficient <- round(nonzerocoefs$coefficient, digits = 4)
    print(nonzerocoefs[order(abs(nonzerocoefs$coefficient), decreasing = TRUE),])
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
#' ("lambda.1se"). 
#' @param type character string. The type of prediction required; the default
#' \code{type = "link"} is on the scale of the linear predictors. Alternatively,
#' for nominal outputs, \code{type = "response"} gives the fitted probabilities
#' and \code{type = "class"} gives the predicted class membership.
#' @param ... currently not used.
#' @details When newdata is not provided, training data included in the specified 
#' object is used.
#' @examples \donttest{
#' airq.ens <- pre(Ozone ~ ., data = airquality[complete.cases(airquality),])
#' predict(airq.ens, newdata = airquality[complete.cases(airquality),])
#' predict(airq.ens)}
#' @import Matrix
#' @export
#' @method predict pre 
predict.pre <- function(object, newdata = NULL, type = "link", 
                         penalty.par.val = "lambda.1se", ...)
{
  if(is.null(newdata)) {
    newdata <- object$modmat
  } else {
    if(!is.data.frame(newdata)) {
      stop("newdata should be a data frame.")
    }
    # check if newdata has the same columns as object$data:
    if(!all(names(object$data) %in% c(names(newdata), object$y_name))) {
      stop("newdata does not contain all predictor variables from the ensemble")
    } else {
      # take all input variables:
      newdata <- newdata[,names(newdata) %in% object$x_names]
      # add temporary y variable to create model.frame:
      newdata[,object$y_name] <- object$data[,object$y_name][1]
      newdata <- model.frame(object$call$formula, newdata)
      # check if all variables have the same levels: 
      if(!all(unlist(sapply(object$data, levels)) == 
              unlist(sapply(newdata[names(object$data)], levels)))) {
        stop("At least one variable in newdata has different levels than the 
             variables used to create the ensemble")
      }
    }
    coefs <- as(coef.glmnet(object$glmnet.fit, s = penalty.par.val), 
                Class = "matrix")
    # if there are rules in the ensemble, they should be evaluated:
    if(object$type != "linear") {
      # get names of rules with nonzero and zero coefficients:
      nonzerorulenames <- names(coefs[coefs!=0,])[grep("rule", names(coefs[coefs!=0,]))]
      zerorulenames <- names(coefs[coefs==0,])[grep("rule", names(coefs[coefs==0,]))]
      if(length(nonzerorulenames) > 0) {
        nonzerorules <- as.character(
          object$rules$description[object$rules$rule %in% nonzerorulenames])
        newrulevars <- data.frame(r1 = as.numeric(with(newdata, eval(parse(
          text = nonzerorules[1])))))
        names(newrulevars) <- nonzerorulenames[1]
        if(length(nonzerorulenames)>1) {
          for(i in 2:length(nonzerorules)) {
            newrulevars[,nonzerorulenames[i]] <- as.numeric(
              with(newdata, eval(parse(text = nonzerorules[i]))))
          }
        }
        # set all rules with zero coefficients to 0:
        if(length(zerorulenames) > 0) {
          for(i in zerorulenames) {
            newrulevars[,i] <- 0
          }
        }
      } else { # only check and assess rules with zero coefficients
        if(length(zerorulenames) > 0) {
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
    if(object$normalize & object$type != "rules") {
      newdata[,names(object$x_scales)] <- scale(
        newdata[,names(object$x_scales)], center = FALSE, scale = object$x_scales)
    }
    if(object$type != "linear") {
      newdata <- cbind(newdata, newrulevars)
    }
    newdata <- MatrixModels::model.Matrix(object$modmat.formula, data = newdata,
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
#' ("lambda.1se"). 
#' @param nvals optional numeric vector of length one. For how many values of x 
#' should the partial dependence plot be created?
#' @param type character string. Type of prediction to be plotted on y-axis. 
#' \code{type = "response"} gives fitted values for continuous outputs and 
#' fitted probabilities for nominal outputs. \code{type = "link"} gives fitted
#' values for continuous outputs and linear predictor values for nominal outputs.
#' @param penalty.par.val character. Penalty parameter criterion to be used for 
#' selecting final model: lambda giving minimum cv error (\code{"lambda.min"}) or 
#' lambda giving cv error that is within 1 standard error of minimum cv error 
#' ("\code{lambda.1se}").
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
#' airq.ens <- pre(Ozone ~ ., data = airquality[complete.cases(airquality),])
#' singleplot(airq.ens, "Temp")}
#' @export
singleplot <- function(object, varname, penalty.par.val = "lambda.1se", 
                       nvals = NULL, type = "response") 
{
  # preliminaries:
  if(length(varname) != 1) {
    stop("A partial dependence plot should be requested for 1 variable")
  }
  if(!is.character(varname)) {
    stop("Specified varname should be of mode character")
  }
  if(is.factor(object$data[,varname]) & !is.null(nvals)) {
    warning("Plot is requested for variable of class factor. Value specified for 
            nvars will be ignored.")
    nvals <- NULL
  }
  
  # Generate expanded dataset:
  if(is.null(nvals)) {
    newx <- unique(object$data[,varname])
  } else {
    newx <- seq(
      min(object$data[,varname]), max(object$data[,varname]), length = nvals)
  }
  exp_dataset <- object$data[rep(row.names(object$data), times = length(newx)),]
  exp_dataset[,varname] <- rep(newx, each = nrow(object$data))

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
#' standard error of minimum cv error ("lambda.1se")?
#' @param phi numeric. See \code{persp()} documentation.
#' @param theta numeric. See \code{persp()} documentation.
#' @param col character. Optional color to be used for surface in 3D plot.
#' @param nvals optional numeric vector of length 2. For how many values of 
#' x1 and x2 should partial dependence be plotted? If \code{NULL}, all observed 
#' values for the two predictor variables specified will be used (see details). 
#' @param ticktype character string. If \code{"simple"} draws just an arrow 
#' parallel to the axis to indicate direction of increase; \code{"detailed"} 
#' draws normal ticks as per 2D plots.
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
#' @examples \donttest{
#' airq.ens <- pre(Ozone ~ ., data = airquality[complete.cases(airquality),])
#' pairplot(airq.ens, c("Temp", "Wind"))}
#' @export
#' @import graphics akima
pairplot <- function(object, varnames, penalty.par.val = "lambda.1se", phi = 45, 
                     theta = 315, col = "cyan", nvals = c(20, 20), ticktype = "detailed",
                     nticks = max(nvals), type = "response", ...) 
{
  # preliminaries:
  if(length(varnames) != 2) {
    stop("Partial dependence should be requested for 2 variables.")
  }
  if(!is.character(varnames)) {
    stop("Specified varname should be of mode character.")
  }
  if(any(sapply(object$data[,varnames], is.factor))) {
    stop("3D partial dependence plots are currently not supported for factors.")
  }
  # generate expanded dataset: 
  if(is.null(nvals)){
    newx1 <- unique(object$data[,varnames[1]])
    newx2 <- unique(object$data[,varnames[2]])
  } else {
    newx1 <- seq(min(object$data[,varnames[1]]), max(object$data[,varnames[1]]), 
                 length = nvals[1])
    newx2 <- seq(min(object$data[,varnames[2]]), max(object$data[,varnames[2]]), 
                 length = nvals[2])
  }
  nobs1 <- length(newx1)
  nobs2 <- length(newx2)
  nobs <- nobs1*nobs2
  exp_dataset <- object$data[rep(row.names(object$data), times = nobs),]
  exp_dataset[,varnames[1]] <- rep(newx1, each = nrow(object$data)*nobs2)
  exp_dataset[,varnames[2]] <- rep(rep(newx2, each = nrow(object$data)), 
                                   times = nobs1)
  
  # get predictions:  
  pred_vals <- predict.pre(object, newdata = exp_dataset, type = type, 
                            penalty.par.val = penalty.par.val)

  # create plot:
  if(is.null(nvals)) nvals <- 3
  xyz <- akima::interp(exp_dataset[,varnames[1]], exp_dataset[,varnames[2]], 
                       pred_vals, duplicate = "mean")
  persp(xyz, xlab = varnames[1], ylab = varnames[2], zlab = "predicted y", 
        phi = phi, theta = theta, col = col, ticktype = ticktype, 
        nticks = nticks, ...)
  
  # to be implemented:
  # type = flag for type of plot when both var1 and var2 are numeric
  # type = "image" => heat map plot
  # type = "persp" => perspective mesh plot
  # type = "contour" => contour plot
  # chgvars = flag for changing plotting relationship when both var1 and var2 are categorical (factors)
  # chgvars = FALSE => plot the partial dependence on the variable (factor) with the most values (levels), for each of the  respective values (levels) of the other variable (factor)
  # chgvars = TRUE => reverse this relationship
  # qntl = trimming factor for plotting numeric variables. Plots are shown for variable values in the range [quantile (qntl) - quantile(1-qntl)]. (Ignored for categorical variables (factors).)
  # nval = maximum number of evaluation points for numeric variables. (Ignored for categorical variables).
  # nav = maximum number of observations used for averaging calculations. (larger values provide higher accuracy with a diminishing return; computation grows linearly with nav)
  # vals1 = vector of names for values (levels) of var1 if it is categorical (factor). (Ignored if var1 is numeric)
  # vals2 = vector of names for values (levels) of var2 if it is categorical (factor). (Ignored if var2 is numeric) 
  # horiz = plot orientation for categorical variable barplots
  # horiz = T/F => do/don't plot bars horizontally
  # las = label orientation flag for categorical variable plots (horiz = F, only)
  # las =1 => horizontal orientation of value (level) names stored in vals1 and/or vals2 (if present).
  # las =2 => vertical orientation of value (level) names stored in vals1 and/or vals2 (if present).
  # cex.names = expansion factor for axis names (bar labels)  for categorical variable barplots 
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
#' @param ... further arguments to be passed to \code{\link[graphics]{barplot}} 
#' (only used when \code{plot = TRUE}).
#' @return A list with two dataframes: $baseimps, giving the importances for 
#' baselearners in the ensemble, and $varimps, giving the importances for 
#' variables that appear and do not appear in the ensemble.
#' @examples \donttest{
#' airq.ens <- pre(Ozone ~ ., data=airquality[complete.cases(airquality),])
#' # calculate global importances:
#' importance(airq.ens) 
#' # calculate local importances (default: over 25% highest predicted values):
#' importance(airq.ens, global = FALSE) 
#' # calculate local importances (custom: over 25% highest predicted values):
#' importance(airq.ens, global = FALSE, quantprobs = c(0, .25))}
#' @export
importance <- function(object, plot = TRUE, ylab = "Importance", 
                       main = "Variable importances", global = TRUE, 
                       quantprobs = c(.75, 1), col = "grey", round = NA, ...) 
{
  ## Step 1: Calculate the importances of the base learners:
  
  # get base learner coefficients (and 'unfactor'):
  coefs <- coef.pre(object, print = FALSE)
  coefs$description <- as.character(coefs$description)
  # give linear terms a description:
  coefs$description[is.na(coefs$description)] <- 
    paste(as.character(coefs$rule)[is.na(coefs$description)], " ", sep = "")
  coefs <- coefs[order(coefs$rule),]
  # Get sds for every baselearner:
  if(global) {
    sds <- c(0, apply(object$modmat, 2, sd, na.rm = TRUE))  
  } else {
    preds <- predict.pre(object, newdata = object$data, type = "response")
    local_modmat <- object$modmat[preds >= quantile(preds, probs = quantprobs[1]) & 
                               preds <= quantile(preds, probs = quantprobs[2]),]
    if(nrow(local_modmat) < 2) {stop("Requested range contains less than 2 
                                observations, importances cannot be calculated")}
    sds <- c(0, apply(local_modmat, 2, sd, na.rm = TRUE)) 
  }
  names(sds)[1] <- "(Intercept)"
  sds <- sds[order(names(sds))]
  if(all(names(sds) != coefs$rule)) {
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
  # For rules, calculate the number of terms in each rule:
  baseimps$nterms <- NA
  for(i in 1:nrow(baseimps)) {
    # If there is no "&" in rule description, there is only 1 term/variable in 
    # the base learner: 
    if(gregexpr("&", baseimps$description)[[i]][1] == -1) {
        baseimps$nterms[i] <- 1 
    } else { # otherwise, the number of terms = the number of &-signs + 1 
      baseimps$nterms[i] <- length(gregexpr("&", baseimps$description)[[i]]) + 1
    }
  }
  # Calculate variable importances:
  varimps <- data.frame(varname = object$x_names, imp = 0)
  # Get importances for rules:
  for(i in 1:nrow(varimps)) { # for every variable:
    # For every baselearner:
    for(j in 1:nrow(baseimps)) {
      # if the variable name appears in the rule:
      if(gregexpr(paste(varimps$varname[i], " ", sep = ""), 
                  baseimps$description[j])[[1]][1] != -1) {
        # then count the number of times it appears in the rule:
        n_occ <- length(
          gregexpr(paste(varimps$varname[i], " ", sep = ""), 
                   baseimps$description[j])[[1]]
        )
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
    if(sum(i == baseimps$modframename, na.rm = TRUE) > 1) {
      varimps$imp[varimps$varname == i] <- sum(varimps$imp[varimps$varname == i],
                      baseimps$imp[i == baseimps$modframename], na.rm = TRUE)
    }
  }
  
  ## Step 3: return (and plot) importances:
  baseimps <- baseimps[baseimps$imp != 0,]
  baseimps <- baseimps[order(baseimps$imp, decreasing = TRUE),]
  varimps <- varimps[order(varimps$imp, decreasing = TRUE),]
  varimps <- varimps[varimps$imp != 0,]
  if(plot == TRUE) {
    barplot(height = varimps$imp, names.arg = varimps$varname, ylab = ylab, 
            main = main, col = col)
  }
  if(!is.na(round)) {
    varimps[,"imp"] <- round(varimps[,"imp"], digits = round)
    baseimps[,c("imp", "coefficient", "sd")] <- round(
      baseimps[,c("imp", "coefficient", "sd")], digits = round)
  }
  return(list(varimps = varimps, baseimps = baseimps[
    baseimps$description != "(Intercept) ", c("description", "imp", "coefficient", "sd")]))
}





#' Compute boostrapped null interaction models
#' 
#' \code{bsnullinteract} calculates null interaction models on bootstrapped 
#' datasets, for deriving a reference distribution of the test statistic 
#' calculated with \code{\link{interact}}.
#' 
#' @param object object of class \code{\link{pre}}.
#' @param nsamp numeric. Number of bootstrapped null interaction models to be 
#' derived.
#' @param seed numeric. Random seed to be used (for reproducability).
#' @return A list of null interaction models, to be used as input for 
#' \code{\link{interact}}.
#' @examples \donttest{
#' airq.ens <- pre(Ozone ~ ., data=airquality[complete.cases(airquality),])
#' nullmods <- bsnullinteract(airq.ens)}
#' @details Computationally intensive. Progress info is printed to command line.
#' @export
bsnullinteract <- function(object, nsamp = 10, seed = 42) {
  # preliminaries:
  set.seed(seed)
  # create call for bootstrapped null model:
  bsnullmodcall <- object$call
  bsnullmodcall$maxdepth <- 1
  bsnullmodcall$verbose <- FALSE  
  # create call for model allowing for interactions, grown on bootstrapped 
  # datasets without interactions:
  bsintmodcall <- bsnullmodcall
  if(is.null(object$call$maxdepth)) {
    bsintmodcall$maxdepth <- 3
  } else {
    bsintmodcall$maxdepth <- object$call$maxdepth
  }
  # compute boostrapped null dataset (i.e., dataset with no interactions):
  bs.ens <- list()
  cat("This may take a while. Computing null model ")
  for(i in 1:nsamp) {
    cat(i, "of", nsamp, "... ")
    # step 1: Take bootstrap sample {x_ip, y_ip}:
    bsdataset <- object$data[sample(1:nrow(object$data), nrow(object$data), 
                                    replace = TRUE),]
    # step 2: Build F_null, a model involving main effects only using {x_ip, y_ip}:
    bsnullmodcall$data <- bsdataset
    bs.ens.null <- eval(bsnullmodcall)
    # step 3: Calculate residuals from predictions y^hat_ip using x_ip and F_null:
    yhatip <- bsdataset[object$y_name] - predict.pre(bs.ens.null, 
                                                      newdata = bsdataset) 
    # step 4: Calculate predictions for original x, using F_null:
    fipx <- predict.pre(bs.ens.null, newdata = object$data)
    # step 5: Calculate ybar, by adding residuals from step 3 to predictions from 
    # step 4:
    bsdataset[,object$y_name] <- yhatip + fipx
    # step 6: Build a model using (x,ybar), using the same procedure as was 
    # originally applied to (x,y):
    bsintmodcall$data <- bsdataset
    bs.ens[[i]] <- eval(match.call(pre, call = bsintmodcall))
  }
  cat("Done!")
  return(bs.ens)
}





# Internal function for calculating H statistic (section 8.1, equation 45):
Hsquaredj <- function(object, varname, k = 10) {
  # Calculate the predicted value F(x) of the full model for each observation:
  preds_x <- predict.pre(object, newdata = object$data)
  # Calculate the expected value of F_j(x_j), over all observed values x_/j,
  # and the expected value of F_/j(x_/j), over all observed values x_j:
  exp_dataset <- object$data[rep(row.names(object$data), 
                                 times = nrow(object$data)),]
  exp_dataset[,varname] <- rep(object$data[,varname], each = nrow(object$data))
  # using predict.pre for a hughe dataset may lead to errors, so split 
  # computations up in 10 parts:
  exp_dataset$ids <- sample(1:k, nrow(exp_dataset), replace = TRUE)
  for(i in 1:k) {  
    cat(".")
    exp_dataset[exp_dataset$ids==i, "yhat"] <- predict.pre(
      object, newdata = exp_dataset[exp_dataset$ids==i,])
  }
  # expected value of F_j(x_j), over all observed values x_/j:
  exp_dataset$i_xj <- rep(1:nrow(object$data), each = nrow(object$data))
  preds_xj <- aggregate(yhat ~ i_xj, data = exp_dataset, FUN = mean)$yhat
  # expected value of F_/j(x_/j), over all observed values x_j:
  exp_dataset$i_xnotj <-  rep(1:nrow(object$data), times = nrow(object$data))
  preds_xnotj <- aggregate(yhat ~ i_xnotj, data = exp_dataset, FUN = mean)$yhat
  # H should be calculated based on centered functions:
  preds_x <- scale(preds_x, center = TRUE, scale = FALSE)
  preds_xj <- scale(preds_xj, center = TRUE, scale = FALSE)
  preds_xnotj <- scale(preds_xnotj, center = TRUE, scale = FALSE)
  if(sum(preds_x^2) > 0) {
    return(sum((preds_x - preds_xj - preds_xnotj)^2) / sum(preds_x^2))
  }
  if(sum(preds_x^2) == 0) {
    warning("The denominator for calculating H squared was equal to zero. It was 
            set to 1e-10 to allow for calculation of H squared")
    return(sum((preds_x - preds_xj - preds_xnotj)^2) / 1e-10)
  }
}





#' Calculate interaction statistics for user-specified variables
#' 
#' \code{interact} calculates a statistic for assessing the presence of 
#' interactions between the input variable(s) specified, and all other input
#' variables.
#' 
#' @param object an object of class \code{\link{pre}}.
#' @param varnames character vector. Names of variables for which interaction 
#' statistics should be calculated.
#' @param k integer. Calculating interaction test statistics is a computationally 
#' intensive, so  calculations are split up in several parts to prevent memory 
#' allocation errors. If a memory allocation error still occurs, increase k.
#' @param nullmods object with bootstrapped null interaction models, resulting 
#' from application of \code{bsnullinteract}.
#' @param plot logical Should the interaction statistics be plotted?
#' @param col character vector of length two. Color for plotting bars used. Only 
#' used when plot = TRUE. Only first element of vector is used if 
#' \code{nullmods = NULL}.
#' @param ylab character string. Label to be used for plotting y-axis.
#' @param ... Additional arguments to be passed to \code{barplot}.
#' @examples 
#' \donttest{
#'  airq.ens <- pre(Ozone ~ ., data=airquality[complete.cases(airquality),])
#'  interact(airq.ens, "Temp")}
#' @details Can be computationally intensive, especially when nullmods is specified.
#' @return If nullmods is not specified, the function returns the interaction 
#' test statistic. If nullmods is specified, the function returns a list, 
#' with elements \code{$H}, which is the test statistic of the interaction 
#' strength, and \code{$nullH}, which is a vector of test statistics of the 
#' interaction in each of the bootstrapped null interaction models. In the barplot,
#' yellow is used for plotting the interaction test statistic. When applicable, 
#' blue is used for the mean in the bootstrapped null models.
#' @export
interact <- function(object, varnames = NULL, nullmods = NULL, k = 10, plot = TRUE, 
                     col = c("yellow", "blue"), ylab = "Interaction strength", 
                     ...)
{
  if(is.null(varnames)) {
    varnames <- as.character(importance(object, plot = FALSE)$varimps$varname)
  }
  if(!all(varnames %in% object$x_names)) {
    stop("Interaction statistics requested for one or more unknown input variables")
  }
  cat("This will take a while (", 
      k*(length(nullmods)+1)*length(varnames), "dots ). ")
  # Calculate H_j for the original dataset:
  H <- rep(NA, length(varnames))
  nullH <- data.frame(matrix(NA, nrow = length(nullmods), ncol = length(varnames)))
  colnames(nullH) <- varnames
  for(i in 1:length(varnames)) {
    H[i] <- Hsquaredj(object = object, varname = varnames[i], k = k)    
    # Calculate mean and sd of H_j for the bootstrapped null models  
    if(!is.null(nullmods)) {
      for(j in 1:length(nullmods)) {
        nullH[j,i] <- Hsquaredj(object = nullmods[[j]], varname = varnames[i], 
                                k = k)  
      }
    }
  }
  cat("\n")
  if(is.null(nullmods)) {
    if(plot) {
      names(H) <- varnames # object$x_names
      barplot(H, col = col[1], ...)
    }
    return(data.frame(varnames = varnames, trainingH2 = H))
  }
  if(!is.null(nullmods)) {
    if(plot) {
      nullmeans <- vector()
      for(i in 1:length(varnames)) {
        nullmeans[i] <- mean(nullH[,i])
      }
      H2s <- matrix(c(H, nullmeans), nrow = 2, byrow = T)
      colnames(H2s) <- varnames
      barplot(H2s, col = col, ylab = ylab, ...)
    }
    return(list(trainingH2 = H, nullH2 = nullH))
  }  
}