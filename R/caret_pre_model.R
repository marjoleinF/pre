#' Model set up for train function of package caret
#' 
#' \code{caret_pre_model} is deprecated and provided for backwards compatibility
#' only. The object provides a model setup for function \code{train} of
#' package caret. It allows for tuning arguments sampfrac, maxdepth, learnrate, 
#' mtry, use.grad and penalty.par.val.
#'  
#' @details 
#' 
#' Object caret_pre_model is deprecated, and only included in package pre for backward 
#' compatibility. Parameters of function \code{pre()} can be tuned by using method 
#' \code{"pre"} in caret's function \code{train()}. See vignette on tuning for more
#' information and examples: \code{vignette("Tuning", package = "pre")} 
#'  
#' @examples 
#' ## Object caret_pre_model is only included in package pre for backward compatibility
#' ## By now, function pre can be optimized in the default way by using the method "pre" 
#' ## in caret's function train(). More information and instructions on tuning parameters
#' ## of function pre() are provided in the vignette about tuning, which can be accessed
#' ## from R by typing:
#' ##
#' ## vignette("Tuning", package = "pre")
#' ##
caret_pre_model <- list(
  library = "pre",
  type = c("Classification", "Regression"),
  parameters = data.frame(parameter = c("sampfrac", "maxdepth", 
                                        "learnrate", "mtry", 
                                        "use.grad", 
                                        "penalty.par.val"),
                          class = c(rep("numeric", times = 4), 
                                    "logical", "character"),
                          label = c("Subsampling Fraction", 
                                    "Max Tree Depth", 
                                    "Shrinkage", 
                                    "# Randomly Selected Predictors",
                                    "Employ Gradient Boosting", 
                                    "Regularization Parameter")),
  grid = function(x, y, len = NULL, search = "grid", 
                  sampfrac = .5, maxdepth = 3L, learnrate = .01, 
                  mtry = Inf, use.grad = TRUE, penalty.par.val = "lambda.1se") {
    if (search == "grid") {
      if (!is.null(len)) {
        maxdepth <- c(3L, 4L, 2L, 5L, 1L, 6:len)[1:len] 
        if (len > 2) {
          sampfrac <- c(.5, .75, 1)
        }
        if (len > 1) {
          penalty.par.val = c("lambda.min", "lambda.1se")
        }
      }
      out <- expand.grid(sampfrac = sampfrac, maxdepth = maxdepth, 
                         learnrate = learnrate, mtry = mtry, 
                         use.grad = use.grad,  
                         penalty.par.val = penalty.par.val)
    } else if (search == "random") {
      out <- data.frame(
        sampfrac = sample(c(.5, .75, 1), size = len, replace = TRUE),
        maxdepth = sample(2L:6L, size = len, replace = TRUE), 
        learnrate = sample(c(0.001, 0.01, 0.1), size = len, replace = TRUE),
        mtry = sample(c(ceiling(sqrt(ncol(x))), ceiling(ncol(x)/3), ncol(x)), size = len, replace = TRUE),
        use.grad = sample(c(TRUE, FALSE), size = len, replace = TRUE),
        penalty.par.val = sample(c("lambda.1se", "lambda.min"), size = len, replace = TRUE))
    }
    return(out)
  },
  fit = function(x, y, wts = NULL, param, lev = NULL, last = NULL, 
                 classProbs, ...) { 
    dat <- if (is.data.frame(x)) x else as.data.frame(x)
    dat$.outcome <- y
    theDots <- list(...)
    if (!any(names(theDots) == "family")) {
      theDots$family <- if (is.factor(y)) {
        if (nlevels(y) == 2L) { 
          "binomial"
        } else {
          "multinomial"
        }
      } else {
        "gaussian"
      }
    }
    
    if(!is.null(wts)) theDots$weights <- wts
    
    modelArgs <- c(
      list(
        formula = as.formula(".outcome ~ ."),
        data = dat,
        sampfrac = param$sampfrac, 
        maxdepth = param$maxdepth, 
        learnrate = param$learnrate, 
        mtry = param$mtry, 
        use.grad = param$use.grad),
        theDots)
    out <- do.call(pre::pre, modelArgs)
    out
  },
  predict = function(modelFit, newdata, submodels = NULL) {
      if (modelFit$family %in% c("gaussian", "mgaussian")) {
        out <- pre:::predict.pre(object = modelFit, 
                                 penalty.par.val = as.character(modelFit$tuneValue$penalty.par.val),
                                 newdata = as.data.frame(newdata))
      } else if (modelFit$family == "poisson") {
        out <- pre:::predict.pre(object = modelFit, 
                                 penalty.par.val = as.character(modelFit$tuneValue$penalty.par.val),
                                 newdata = as.data.frame(newdata), type = "response")
      } else {
        out <- factor(pre:::predict.pre(object = modelFit, 
                                        penalty.par.val = as.character(modelFit$tuneValue$penalty.par.val),
                                        newdata = as.data.frame(newdata), type = "class"))      
      }
    if (!is.null(submodels)) {
      tmp <- list()
      for (i in seq(along.with = submodels$penalty.par.val)) {
        if (modelFit$family %in% c("gaussian", "mgaussian")) {
          tmp[[i]] <- pre:::predict.pre(object = modelFit, 
                                        newdata = as.data.frame(newdata), 
                                        penalty.par.val = as.character(submodels$penalty.par.val[i])) 
        } else if (modelFit$family == "poisson") {
          tmp[[i]] <- pre:::predict.pre(object = modelFit, 
                                        newdata = as.data.frame(newdata), 
                                        type = "response",
                                        penalty.par.val = as.character(submodels$penalty.par.val[i]))
        } else {
          tmp[[i]] <- factor(pre:::predict.pre(object = modelFit, 
                                               newdata = as.data.frame(newdata), 
                                               type = "class",
                                               penalty.par.val = as.character(submodels$penalty.par.val[i])))      
        }
      }
      out <- c(list(out), tmp)
    }
    out
  },
  prob = function(modelFit, newdata, submodels = NULL) {
      probs <- pre:::predict.pre(object = modelFit, 
                                 penalty.par.val = as.character(modelFit$tuneValue$penalty.par.val),
                                 newdata = as.data.frame(newdata), 
                                 type = "response")
      if (is.null(ncol(probs)) || ncol(probs) == 1) {
        probs <- data.frame(1 - probs, probs)
        colnames(probs) <- levels(modelFit$data[,modelFit$y_names])
      }
    if (!is.null(submodels)) {
      tmp <- list()
      for (i in seq(along.with = submodels$penalty.par.val)) {
        tmp[[i]] <- pre:::predict.pre(object = modelFit, 
                                      newdata = as.data.frame(newdata), 
                                      type = "response",
                                      penalty.par.val = as.character(submodels$penalty.par.val[i]))
        if (is.null(ncol(tmp[[i]])) || ncol(tmp[[i]]) == 1) {
          tmp[[i]] <- data.frame(1 - tmp[[i]], tmp[[i]])
          colnames(tmp[[i]]) <- levels(modelFit$data[ , modelFit$y_names])
        }
      }
      probs <- c(list(probs), tmp)
    }
    probs
  },
  sort = function(x) {
    ordering <- order(x$maxdepth, # lower values are simpler
                      x$use.grad, # TRUE employs ctree (vs ctree), so simplest
                      max(x$mtry) - x$mtry, # higher values yield more similar tree, so simpler
                      x$sampfrac != 1L, # subsampling yields simpler trees than bootstrap sampling
                      x$learnrate, # lower learnrates yield more similar trees, so simpler
                      decreasing = FALSE)
    x[ordering,]
  },
  loop = function(fullGrid) {
    
    loop_rows <- rownames(unique(fullGrid[,-which(names(fullGrid) == "penalty.par.val")]))
    loop <- fullGrid[rownames(fullGrid) %in% loop_rows, ]
    
    submodels <- list()

    for (i in 1:nrow(loop)) {
      lambda_vals <- character()

      for (j in 1:nrow(fullGrid)) {
        if (all(loop[i, -which(colnames(loop) == "penalty.par.val")] ==
                fullGrid[j, -which(colnames(fullGrid) == "penalty.par.val")])) {
          lambda_vals <- c(lambda_vals, as.character(fullGrid[j, "penalty.par.val"]))
        }
      }
      lambda_vals <- lambda_vals[-which(lambda_vals == loop$penalty.par.val[i])]
      submodels[[i]] <- data.frame(penalty.par.val = lambda_vals)
    }
    list(loop = loop, submodels = submodels)
  },
  levels = NULL,
  tag = c("Rule-Based Model", "Tree-Based Model", "L1 regularization", "Bagging", "Boosting"),
  label = "Prediction Rule Ensembles",
  predictors = NULL,
  varImp = NULL, 
  oob = NULL,
  notes = NULL,
  check = NULL
)