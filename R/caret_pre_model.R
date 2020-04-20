#' Model set up for train function of package caret
#' 
#' \code{caret_pre_model} provides a model setup for function \code{train} of
#' package caret. It allows for tuning arguments sampfrac, maxdepth, learnrate, 
#' mtry, use.grad and penalty.par.val.
#'  
#' @details 
#' 
#' When tuning parameters of \code{pre()} with caret's \code{train()}
#' function, always use the default S3 method (i.e., specify predictors and response
#' variables through arguments \code{x} and \code{y}). When \code{train.formula()}, 
#' is used (i.e., if \code{formula} and \code{data} arguments are specified),
#' \code{train} will internally call \code{model.matrix()} on \code{data}, which 
#' will code all categorical (factor) predictor variables as dummy variables,
#' which will yield a different result than inputting the original factors, for most 
#' tree-based methods.
#' 
#' \code{caret_pre__model$parameters} provides an overview of the parameters that
#' can be tuned for function \code{pre} using \code{caret}. \code{caret_pre_model$grid}
#' provides a function for creating a tuning grid (see Examples below).
#' 
#' @examples \dontrun{
#'  
#' library("caret")
#'
#' ## Prepare data:
#' airq <- airquality[complete.cases(airquality),]
#' y <- airq$Ozone
#' x <- airq[,-1]
#' 
#' ## Apply caret with only pre's default settings (trControl and ntrees argument
#' ## are employed here only to reduce computation time):
#' 
#' set.seed(42)
#' prefit1 <- train(x = x, y = y, method = caret_pre_model,
#'                  trControl = trainControl(number = 1),
#'                  ntrees = 25L)
#' prefit1
#' 
#' ## Create custom tuneGrid:
#' set.seed(42)
#' tuneGrid <- caret_pre_model$grid(x = x, y = y,
#'                                  maxdepth = 3L:5L,
#'                                  learnrate = c(.01, .1),
#'                                  penalty.par.val = c("lambda.1se", "lambda.min"))
#' tuneGrid
#' ## Apply caret (again, ntrees and trControl set only to reduce computation time):
#' prefit2 <- train(x = x, y = y, method = caret_pre_model,
#'                  trControl = trainControl(number = 1),
#'                  tuneGrid = tuneGrid, ntrees = 25L)
#' prefit2
#' 
#' ## Get best tuning parameter values:
#' prefit2$bestTune
#' ## Get predictions from model with best tuning parameters:
#' predict(prefit2, newdata = x[1:10, ])
#' plot(prefit2)
#'
#' ## Obtain tuning grid through random search over the tuning parameter space:
#' set.seed(42)
#' tuneGrid2 <- caret_pre_model$grid(x = x, y = y, search = "random", len = 10)
#' tuneGrid2
#' set.seed(42)
#' prefit3 <- train(x = x, y = y, method = caret_pre_model,
#'                  trControl = trainControl(number = 1, verboseIter = TRUE),
#'                  tuneGrid = tuneGrid2, ntrees = 25L)
#' prefit3
#' 
#' ## Count response:
#' set.seed(42)
#' prefit4 <- train(x = x, y = y, method = caret_pre_model,
#'                  trControl = trainControl(number = 1),
#'                  ntrees = 25L, family = "poisson")
#' prefit4
#' 
#' ## Binary factor response:
#' y_bin <- factor(airq$Ozone > mean(airq$Ozone))
#' set.seed(42)
#' prefit5 <- train(x = x, y = y_bin, method = caret_pre_model,
#'                  trControl = trainControl(number = 1),
#'                  ntrees = 25L, family = "binomial")
#' prefit5
#' 
#' ## Factor response with > 2 levels:
#' x_multin <- airq[,-5]
#' y_multin <- factor(airq$Month)
#' set.seed(42)
#' prefit6 <- train(x = x_multin, y = y_multin, method = caret_pre_model,
#'                  trControl = trainControl(number = 1),
#'                  ntrees = 25L, family = "multinomial")
#' prefit6
#' }
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

    #pre(formula = as.formula(".outcome ~ ."), data = dat, weights = weights, 
    #    sampfrac = param$sampfrac, maxdepth = param$maxdepth, 
    #    learnrate = param$learnrate, mtry = param$mtry, 
    #    use.grad = param$use.grad, ...)
  },
  predict = function(modelFit, newdata, submodels = NULL) {
    #if (is.null(submodels)) {
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
    #} else {
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
    #if (is.null(submodels)) {
      probs <- pre:::predict.pre(object = modelFit, 
                                 penalty.par.val = as.character(modelFit$tuneValue$penalty.par.val),
                                 newdata = as.data.frame(newdata), 
                                 type = "response")
      # For binary classification, create matrix:    
      if (is.null(ncol(probs)) || ncol(probs) == 1) {
        probs <- data.frame(1 - probs, probs)
        colnames(probs) <- levels(modelFit$data[,modelFit$y_names])
      }
    #} else {
    if (!is.null(submodels)) {
      tmp <- list()
      for (i in seq(along.with = submodels$penalty.par.val)) {
        tmp[[i]] <- pre:::predict.pre(object = modelFit, 
                                      newdata = as.data.frame(newdata), 
                                      type = "response",
                                      penalty.par.val = as.character(submodels$penalty.par.val[i]))
        # For binary classification, create matrix:    
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
    
    # loop should provide a grid containing models that can
    # be looped over for tuning penalty.par.val
    loop_rows <- rownames(unique(fullGrid[,-which(names(fullGrid) == "penalty.par.val")]))
    loop <- fullGrid[rownames(fullGrid) %in% loop_rows, ]
    
    ## submodels should be a list and length(submodels == nrow(loop)
    ## each element of submodels should be a data.frame with column penalty.par.val, with a row for every value to loop over
    submodels <- list()
    ## for every row of loop:
    for (i in 1:nrow(loop)) {
      lambda_vals <- character()
      ## check which rows in fullGrid without $penalty.par.val are equal to
      ## rows in loop without $penalty.par.val
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
  levels = NULL, #function(x) { levels(x$data[,x$y_names]) },
  tag = c("Rule-Based Model", "Tree-Based Model", "L1 regularization", "Bagging", "Boosting"),
  label = "Prediction Rule Ensembles",
  predictors = NULL, #function(x, ...) { 
    #if (x$family %in% c("gaussian", "poisson", "binomial")) {
    #  return(suppressWarnings(importance(x, plot = FALSE, 
    #                                     penalty.par.val = object$tuneValue$penalty.par.val,
    #                                     ...)$varimps$varname))
    #} else {
    #  warning("Reporting the predictors in the model is not yet available for multinomial and multivariate responses")
    #  return(NULL)
    #}
  #},
  varImp = NULL, #function(x, ...) {
    #if (x$family %in% c("gaussian","binomial","poisson")) {
    #  varImps <- pre:::importance(x, plot = FALSE, 
    #                              penalty.par.val = object$tuneValue$penalty.par.val, 
    #                              ...)$varimps
    #  varnames <- varImps$varname
    #  varImps <- data.frame(Overall = varImps$imp)
    #  rownames(varImps) <- varnames  
    #  return(varImps)
    #} else {
    #  warning("Variable importances cannot be calculated for multinomial or mgaussian family")
    #  return(NULL)
    #}
  #},
  oob = NULL,
  notes = NULL,
  check = NULL
)


