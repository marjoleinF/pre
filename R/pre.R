utils::globalVariables("%dopar%")

#' Derive a prediction rule ensemble
#'
#' \code{pre} derives a sparse ensemble of rules and/or linear functions for 
#' prediction of a continuous, binary, count, multinomial, multivariate 
#' continuous or survival response.
#' 
#' @param formula a symbolic description of the model to be fit of the form 
#' \code{y ~ x1 + x2 + ...+ xn}. Response (left-hand side of the formula) 
#' should be of class numeric (for \code{family = "gaussian"} or
#' \code{"mgaussian"}), integer (for \code{family = "poisson"}), factor (for 
#' \code{family = "binomial"} or \code{"multinomial"}). See Examples below. 
#' Note that the minus sign (\code{-}) may not be used in the formula to omit
#' the intercept or variables in \code{data}, and neither should  \code{+ 0} 
#' be used to omit the intercept. To omit the intercept from the final ensemble, 
#' add \code{intercept = FALSE} to the call. To omit variables from the final
#' ensemble, make sure they are excluded from \code{data}.
#' @param data \code{data.frame} containing the variables in the model. Response 
#' must be of class \code{factor} for classification, \code{numeric} for (count) 
#' regression, \code{Surv} for survival regression. Input variables must be of 
#' class numeric, factor or ordered factor. Otherwise, \code{pre} will attempt
#' to recode.
#' @param family specifies a glm family object. Can be a character string (i.e., 
#' \code{"gaussian"}, \code{"binomial"}, \code{"poisson"}, \code{"multinomial"}, 
#' \code{"cox"} or \code{"mgaussian"}), or a corresponding family object 
#' (e.g., \code{gaussian}, \code{binomial} or \code{poisson}, see 
#' \code{\link[stats]{family}}). Specification of argument \code{family} is 
#' strongly advised but not required. If \code{family} is not specified, 
#' Otherwise, the program will try to make an informed guess, based on the 
#' class of the response variable specified in \code{formula}. als see Examples 
#' below. 
#' @param use.grad logical. Should gradient boosting with regression trees be
#' employed when \code{learnrate > 0}? If \code{TRUE}, use trees fitted by 
#' \code{\link[partykit]{ctree}} or \code{\link[rpart]{rpart}} as in Friedman 
#' (2001), but without the line search. If \code{use.grad = FALSE}, 
#' \code{\link[partykit]{glmtree}} instead of \code{\link[partykit]{ctree}} 
#' will be employed for rule induction, yielding longer computation times, 
#' higher complexity, but possibly higher predictive accuracy. See Details for 
#' supported combinations of \code{family}, \code{use.grad} and \code{learnrate}.
#' @param weights optional vector of observation weights to be used for 
#' deriving the ensemble.
#' @param type character. Specifies type of base learners to include in the 
#' ensemble. Defaults to \code{"both"} (initial ensemble will include both rules 
#' and linear functions). Other option are \code{"rules"} (prediction 
#' rules only) or \code{"linear"} (linear functions only).
#' @param sampfrac numeric value \eqn{> 0} and \eqn{\le 1}. Specifies  
#' the fraction of randomly selected training observations used to produce each 
#' tree. Values \eqn{< 1} will result in sampling without replacement (i.e., 
#' subsampling), a value of 1 will result in sampling with replacement 
#' (i.e., bootstrap sampling). Alternatively, a sampling function may be supplied, 
#' which should take arguments \code{n} (sample size) and \code{weights}. 
#' @param maxdepth positive integer. Maximum number of conditions in rules. 
#' If \code{length(maxdepth) == 1}, it specifies the maximum depth of 
#' of each tree grown. If \code{length(maxdepth) == ntrees}, it specifies the
#' maximum depth of every consecutive tree grown. Alternatively, a random
#' sampling function may be supplied, which takes argument \code{ntrees} and 
#' returns integer values. See also \code{\link{maxdepth_sampler}}.
#' @param learnrate numeric value \eqn{> 0}. Learning rate or boosting parameter.
#' @param mtry positive integer. Number of randomly selected predictor variables for 
#' creating each split in each tree. Ignored when \code{tree.unbiased=FALSE}.
#' @param ntrees positive integer value. Number of trees to generate for the 
#' initial ensemble.
#' @param confirmatory character vector. Specifies one or more confirmatory terms 
#' to be included in the final ensemble. Linear terms can be specified as the 
#' name of a predictor variable included in \code{data}, rules can be specified
#' as, for example, \code{"x1 > 6 & x2 <= 8"}, where x1 and x2 should be names
#' of variables in \code{data}. Terms thus specified will be included in the
#' final ensemble, as their coefficients will not be penalized in the estimation.
#' @param removeduplicates logical. Remove rules from the ensemble which are 
#' identical to an earlier rule?
#' @param removecomplements logical. Remove rules from the ensemble which are
#' identical to (1 - an earlier rule)? 
#' @param winsfrac numeric value \eqn{> 0} and \eqn{\le 0.5}. Quantiles of data 
#' distribution to be used for 
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
#' @param ordinal logical. Should ordinal variables (i.e., ordered factors) be
#' treated as continuous for generating rules? \code{TRUE} (the default)
#' generally yields simpler rules, shorter computation times and better 
#' generalizability of the final ensemble. 
#' @param nfolds positive integer. Number of cross-validation folds to be used for 
#' selecting the optimal value of the penalty parameter \eqn{\lambda} in selecting
#' the final ensemble.
#' @param verbose logical. Should progress be printed to the command line?
#' @param par.init logical. Should parallel \code{foreach} be used to generate 
#' initial ensemble? Only used when \code{learnrate == 0}. Note: Must register 
#' parallel beforehand, such as doMC or others. Furthermore, setting 
#' \code{par.init = TRUE} will likely only increase computation time for smaller 
#' datasets.
#' @param par.final logical. Should parallel \code{foreach} be used to perform cross 
#' validation for selecting the final ensemble? Must register parallel beforehand, 
#' such as doMC or others.
#' @param tree.control list with control parameters to be passed to the tree 
#' fitting function, generated using \code{\link[partykit]{ctree_control}},
#' \code{\link[partykit]{mob_control}} (if \code{use.grad = FALSE}), or 
#' \code{\link[rpart]{rpart.control}} (if \code{tree.unbiased = FALSE}).
#' @param tree.unbiased logical. Should an unbiased tree generation algorithm 
#' be employed for rule generation? Defaults to \code{TRUE}, if set to 
#' \code{FALSE}, rules will be generated employing the CART algorithm
#' (which suffers from biased variable selection) as implemented in 
#' \code{\link[rpart]{rpart}}. See details below for possible combinations 
#' with \code{family}, \code{use.grad} and \code{learnrate}.
#' @param sparse logical. Should sparse design matrices be used? Likely improves
#' computation times for large datasets.
#' @param ... Additional arguments to be passed to
#' \code{\link[glmnet]{cv.glmnet}}.
#' 
#' @details Note that obervations with missing values will be removed prior to 
#' analysis.
#' 
#' In some cases, duplicated variable names may appear in the model.
#' For example, the first variable is a factor named 'V1' and there are also
#' variables named 'V10' and/or 'V11' and/or 'V12' (etc). Then for 
#` selecting the final ensemble, if linear terms are also included,
#' for the binary factor V1, dummy contrast variables will be created, named 
#' 'V10', 'V11', 'V12' (etc). As should be clear from this example, this yields 
#' duplicated variable names, which may yield problems, for example in the 
#' calculation of predictions and importances, later on. This can be prevented 
#' by renaming factor variables with numbers in their name, prior to analysis.
#' 
#' The table below provides an overview of combinations of response 
#' variable types, \code{use.grad}, \code{tree.unbiased} and
#' \code{learnrate} settings that are supported, and the tree induction 
#' algorithm that will be employed as a result:
#' 
#' \tabular{lccccc}{
#' \strong{use.grad} \tab \strong{tree.unbiased} \tab \strong{learnrate} \tab \strong{family} \tab \strong{tree alg.} \tab \strong{Response variable format} \cr
#' \cr
#' TRUE	\tab TRUE	\tab 0 \tab gaussian	  \tab ctree\tab Single, numeric (non-integer) \cr
#' TRUE	\tab TRUE	\tab 0 \tab mgaussian	  \tab ctree\tab Multiple, numeric (non-integer) \cr
#' TRUE	\tab TRUE	\tab 0 \tab binomial	  \tab ctree\tab Single, factor with 2 levels \cr
#' TRUE	\tab TRUE	\tab 0 \tab multinomial	\tab ctree\tab Single, factor with >2 levels \cr
#' TRUE	\tab TRUE	\tab 0 \tab poisson	    \tab ctree\tab Single, integer \cr
#' TRUE \tab TRUE \tab 0 \tab cox         \tab ctree\tab Object of class 'Surv' \cr
#' \cr
#' TRUE	\tab TRUE	\tab >0 \tab 	gaussian	  \tab ctree \tab Single, numeric (non-integer) \cr
#' TRUE	\tab TRUE	\tab >0	\tab mgaussian	  \tab ctree \tab Multiple, numeric (non-integer) \cr
#' TRUE	\tab TRUE	\tab >0	\tab binomial	  \tab ctree  \tab Single, factor with 2 levels \cr
#' TRUE	\tab TRUE	\tab >0	\tab multinomial	\tab ctree \tab Single, factor with >2 levels \cr
#' TRUE	\tab TRUE	\tab >0	\tab poisson	    \tab ctree  \tab Single, integer \cr
#' TRUE \tab TRUE \tab >0 \tab cox         \tab ctree\tab Object of class 'Surv' \cr
#' \cr
#' FALSE \tab TRUE \tab 0 \tab gaussian	  \tab glmtree \tab Single, numeric (non-integer) \cr
#' FALSE \tab TRUE \tab 0 \tab binomial	  \tab glmtree \tab Single, factor with 2 levels \cr
#' FALSE \tab TRUE \tab 0 \tab poisson	    \tab glmtree \tab Single, integer \cr
#' \cr
#' FALSE \tab TRUE \tab >0 \tab gaussian	  \tab glmtree \tab Single, numeric (non-integer) \cr
#' FALSE \tab TRUE \tab >0 \tab binomial	  \tab glmtree \tab Single, factor with 2 levels \cr
#' FALSE \tab TRUE \tab >0 \tab poisson	    \tab glmtree \tab Single, integer \cr
#' \cr
#' TRUE	\tab FALSE \tab 0 \tab gaussian	  \tab rpart \tab Single, numeric (non-integer) \cr
#' TRUE	\tab FALSE \tab 0 \tab binomial	  \tab rpart \tab Single, factor with 2 levels \cr
#' TRUE	\tab FALSE \tab 0 \tab multinomial	\tab rpart \tab Single, factor with >2 levels \cr
#' TRUE	\tab FALSE \tab 0 \tab poisson	    \tab rpart \tab Single, integer \cr
#' TRUE \tab FALSE \tab 0 \tab cox         \tab rpart\tab Object of class 'Surv' \cr
#' \cr
#' TRUE \tab FALSE	\tab >0 \tab gaussian	  \tab rpart \tab Single, numeric (non-integer) \cr
#' TRUE \tab FALSE	\tab >0 \tab binomial	  \tab rpart \tab Single, factor with 2 levels \cr
#' TRUE \tab FALSE	\tab >0 \tab poisson	  \tab rpart \tab Single, integer \cr
#' TRUE \tab FALSE \tab >0 \tab cox         \tab rpart \tab Object of class 'Surv'
#' }
#' 
#' @note Parts of the code for deriving rules from the nodes of trees was copied 
#' with permission from an internal function of the \code{partykit} package, written
#' by Achim Zeileis and Torsten Hothorn.
#' 
#' @return An object of class \code{pre}. It contains the initial ensemble of 
#' rules and/or linear terms and a range of possible final ensembles. 
#' By default, the final ensemble employed by all other
#' methods and functions in package \code{pre} is selected using the 'minimum
#' cross validated error plus 1 standard error' criterion. All functions and 
#' methods for objects of class \code{pre} take a \code{penalty.parameter.value} 
#' argument, which can be used to select a different criterion.
#' 
#' @examples \donttest{## Fit pre to a continuous response:
#' airq <- airquality[complete.cases(airquality), ]
#' set.seed(42)
#' airq.ens <- pre(Ozone ~ ., data = airq)
#' airq.ens
#' 
#' ## Fit pre to a binary response:
#' airq2 <- airquality[complete.cases(airquality), ]
#' airq2$Ozone <- factor(airq2$Ozone > median(airq2$Ozone))
#' set.seed(42)
#' airq.ens2 <- pre(Ozone ~ ., data = airq2, family = "binomial")
#' airq.ens2
#' 
#' ## Fit pre to a multivariate continuous response:
#' airq3 <- airquality[complete.cases(airquality), ] 
#' set.seed(42)
#' airq.ens3 <- pre(Ozone + Wind ~ ., data = airq3, family = "mgaussian")
#' airq.ens3
#' 
#' ## Fit pre to a multinomial response:
#' set.seed(42)
#' iris.pre <- pre(Species ~ ., data = iris, family = "multinomial")
#' iris.pre
#' 
#' ## Fit pre to a survival response:
#' library("survival")
#' lung <- lung[complete.cases(lung), ]
#' set.seed(42)
#' lung.ens <- pre(Surv(time, status) ~ . - sex, data = lung, family = "cox")
#' lung.ens
#' 
#' ## Fit pre to a count response:
#' ## Generate random data (partly based on Dobson (1990) Page 93: Randomized 
#' ## Controlled Trial):
#' counts <- rep(as.integer(c(18, 17, 15, 20, 10, 20, 25, 13, 12)), times = 10)
#' outcome <- rep(gl(3, 1, 9), times = 10)
#' treatment <- rep(gl(3, 3), times = 10)
#' noise1 <- 1:90
#' set.seed(1)
#' noise2 <- rnorm(90)
#' countdata <- data.frame(treatment, outcome, counts, noise1, noise2)
#' set.seed(42)
#' count.ens <- pre(counts ~ ., data = countdata, family = "poisson")
#' count.ens}
#' @import glmnet partykit datasets
#' @seealso \code{\link{print.pre}}, \code{\link{plot.pre}}, 
#' \code{\link{coef.pre}}, \code{\link{importance}}, \code{\link{predict.pre}}, 
#' \code{\link{interact}}, \code{\link{cvpre}} 
#' @references
#' 
#' Fokkema, M. (2018). Fitting prediction rule ensembles with R package pre.
#' \url{https://arxiv.org/abs/1707.07149}.
#' 
#' Friedman, J. H. (2001). Greedy function approximation: a gradient boosting 
#' machine. \emph{The Annals of Applied Statistics, 29}(5), 1189-1232.
#' 
#' Friedman, J. H., & Popescu, B. E. (2008). Predictive learning via rule 
#' ensembles. \emph{The Annals of Applied Statistics, 2}(3), 916-954.
#' 
#' Hothorn, T., & Zeileis, A. (2015). partykit: A modular toolkit for recursive 
#' partytioning in R. \emph{Journal of Machine Learning Research, 16}, 3905-3909.
#' 
#' @export
pre <- function(formula, data, family = gaussian,
                use.grad = TRUE, weights, type = "both", sampfrac = .5, 
                maxdepth = 3L, learnrate = .01, mtry = Inf, ntrees = 500,
                confirmatory = NULL,
                removecomplements = TRUE, removeduplicates = TRUE, 
                winsfrac = .025, normalize = TRUE, standardize = FALSE,
                ordinal = TRUE, nfolds = 10L, tree.control, tree.unbiased = TRUE, 
                verbose = FALSE, par.init = FALSE, par.final = FALSE, 
                sparse = FALSE, ...) { 
  
  
  #####################
  ## Check arguments ##
  #####################
  
  ## Save call:
  cl <- match.call()
  
  ## Check if proper formula argument is specified 
  ## and check if glmertree should be employed:
  if (!(inherits(formula, "formula"))) {
    stop("Argument 'formula' should specify and object of class 'formula'.\n")
  } else {
    if (length(as.Formula(formula))[[2]] == 3L) { # then right-hand side of regression formula consists of three parts and glmertree should be employed
      formula <- as.Formula(formula)
      use_glmertree <- TRUE
      if (formula[[3]][[2]][[2]] == 1) { # check if intercept is specified as regressor for linear model
          if (use.grad || learnrate > 0) {
            stop("When specifying a formula with three-part right-hand side, argument 'use.grad' should be set to FALSE and 'learnrate' to 0.\n", immediate. = TRUE)
          }
        } else {
          stop("When specifying a three-part right-hand side, the first part of the right-hand side should consist of an intercept only, e.g., y ~ 1 | cluster | x1 + x2 + x3.\n")
        }
      } else {
      use_glmertree <- FALSE
      }
    ## Check if dot and functions are simultaneously used in formula
    form <- as.character(formula[3])
    for (i in names(data)) {
      form <- gsub(pattern = i, replacement = "", x = form)
    }
    if (any(grepl(".", form, fixed = TRUE))) {
      if (any(grepl("(", form, fixed = TRUE))) {
        if (any(grepl(")", form, fixed = TRUE))) {
          warning("Argument 'formula' contains both one or more functions of predictor variables, as well as a dot ('.'), which should be avoided. Model fitting may fail, and/or both the original variable(s) and their functions may be included as predictor variables.\n", immediate. = TRUE)  
        }
      }
    }
    if (any(grepl("-", as.character(formula), fixed = TRUE))) {
      warning("Argument 'formula' contains a minus sign. Note that the minus sign should not be used to omit the intercept or variables from the ensemble.\n", immediate. = TRUE)
    }
  }

  ## Check if proper data argument is specified:
  if (!is.data.frame(data)) {
    stop("Argument 'data' should specify a data frame.\n")
  }

  ## Check and set up family argument: 
  if (length(family) > 1L) {
    warning("Argument 'family' has length > 1, only first element will be used.\n")
    family <- family[1L]
  }
  if (is.function(family)) { 
    family <- family()
    if (inherits(family, "family")) {
      link <- family$link
      family <- family$family
      if (family == "gaussian" && link != "identity") {
        warning("The link function specified is currently not supported; identity link will be employed.\n", immediate. = TRUE)
      } else if (family == "binomial" && link != "logit") {
        warning("The link function specified is currently not supported; logit link will be employed.\n", immediate. = TRUE)
      } else if (family == "poisson" && link != "log") {
        warning("The link function specified is currently not supported; log link will be employed.\n", immediate. = TRUE)
      }
    } else {
      stop("Argument 'family' should specify a family object.\n")
    }
  }
  if (is.character(family)) {
    if (!(family %in% c("gaussian", "binomial", "poisson", "mgaussian", "multinomial", "cox"))) {
      stop("Argument 'family' should be equal to 'gaussian', 'binomial', 'poisson', 'multinomial', 'mgaussian', 'cox', or a corresponding family object.\n")
    }
  } else {
    stop("Argument 'family' should be equal to 'gaussian', 'binomial', 'poisson', 'multinomial', 'mgaussian', 'cox', or a corresponding family object.\n")
  }
  
  ## Check if proper weights argument is specified, if specified:
  if (missing(weights)) {
    weights <- rep(1L, times = nrow(data))
  } else if (length(weights) != nrow(data)) {
      warning("Length of argument 'weights' is not equal to nrow(data).\n", immediate. = TRUE)
  }
  
  ## Check if proper type argument is specified:
  if (!(length(type) == 1L && type %in% c("rules", "both", "linear"))) {
    stop("Argument 'type' should be 'rules', 'linear' or 'both'.\n")
  }
  
  ## Check if proper sampfrac argument is specified:
  if (!is.function(sampfrac)) {
    if (!(length(sampfrac) == 1 && is.numeric(sampfrac) && sampfrac > 0 && 
          sampfrac <= 1)) {
      stop("Argument 'sampfrac' should be a single numeric value > 0 and <= 1, or a sampling function.\n")
    }
  }

  ## Check if proper maxdepth argument is specified:
  if (is.function(maxdepth)) {
    maxdepth <- maxdepth(ntrees = ntrees)
  } else if (!is.numeric(maxdepth)) {
    stop("Argument 'maxdepth' should be either a numeric vector of length 1 or ntrees, or a random number generating function.\n")
  } else if (!(length(maxdepth) %in% c(1L, ntrees))) {
    warning("Argument 'maxdepth' should be either a numeric vector of length 1 or ntrees, only first value of maxdepth will be used.\n")
    maxdepth <- maxdepth[1L]
  } 
  if (!all(maxdepth > 0)) {
    stop("All values of maxdepth should be > 0.\n")
  } 
  if (!all(maxdepth == suppressWarnings(as.integer(maxdepth)) | is.infinite(maxdepth))) {
    stop("Argument 'maxdepth' should consist of  of integer values or Inf), or a random number generating function.\n")
  }
  
  ## Check if proper learnrate argument is specified:
  if (!(length(learnrate) == 1L && is.numeric(learnrate) && 
        (learnrate >= 0 || learnrate <= 1))) {
    stop("Argument 'learnrate' shoud be a single numeric value >= 0 and <= 1.\n")
  }
  
  ## Check if proper mtry argument is specified:
  if (!(length(mtry) == 1L && mtry > 0 && 
        (mtry == suppressWarnings(as.integer(mtry)) || is.infinite(mtry)))) {
    stop("Argument 'mtry' should be a single integer value, or Inf.\n")
  }
  
  ## Check if proper ntrees argument is specified:
  if (!(length(ntrees) == 1L && ntrees == as.integer(ntrees) && ntrees > 0)) {
    stop("Argument 'ntrees' should be a single positive integer.\n")
  }
  
  ## Check if proper winsfrac argument is specified:
  if (!(length(winsfrac == 1L) && is.numeric(winsfrac) && winsfrac >= 0 && 
        winsfrac < 1)) {
    stop("Argument 'winsfrac' should be a numeric value >= 0 and < 1.\n")
  }

  ## Check if proper nfolds argument is specified:
  if (!(length(nfolds) == 1L && nfolds > 0 && nfolds == as.integer(nfolds))) {
    stop("Argument 'nfolds' should be a positive integer.\n")
  }
  
  ## Check if proper confirmatory argument is specified:
  if (!is.null(confirmatory)) {
    if (!is.character(confirmatory)) {
      stop("Argument 'confirmatory' should specify a character vector.\n")
    }
  }
  
  ## Check if logical arguments of length 1 are properly specified:
  is_logical_and_length_one <- function(x) {is.logical(x) && length(x) == 1L}
  for (i in c(use.grad, removeduplicates, removecomplements, normalize, 
              standardize, ordinal, verbose, tree.unbiased, par.init, 
              par.final)) {
    if (!is_logical_and_length_one(i)) {
      stop("Argument ", i, "should be TRUE of FALSE.\n")
    }
  }
  
  if (par.final || par.init) {
    if(!requireNamespace("foreach")) {
      warning("Parallel computation requires package foreach. Arguments 'par.init' and 'par.final' are set to FALSE.\n")   
      par.init <- par.final <- FALSE
    }
  }

  ## Check if proper tree.control argument is specified:
  if (missing(tree.control)) {
    if (tree.unbiased && (use.grad || !use_glmertree)) {
      tree.control <- ctree_control(maxdepth = maxdepth[1], mtry = mtry)
    } else if (tree.unbiased && !use.grad) {
      tree.control <- mob_control(maxdepth = maxdepth[1] + 1, mtry = mtry)
    } else if (!tree.unbiased) {
      if (any(maxdepth > 29)) {
        maxdepth[maxdepth > 29] <- 29L
        warning("If tree.unbiased = FALSE, max(maxdepth) is 29.\n")
      }
      tree.control <- rpart.control(maxdepth = maxdepth[1L])
      if (!is.infinite(mtry)) {
        warning("Value specified for mtry will be ignored if tree.unbiased = FALSE.\n")
      }
    }
  } else {
    if (!is.list(tree.control)) {
      stop("Argument 'tree.control' should be a list of control parameters.\n")
    }
    if (use.grad && tree.unbiased && !use_glmertree) {
      if (!setequal(names(ctree_control()), names(tree.control))) {
        stop("Argument 'tree.control' should be a list containing named elements ", paste(names(ctree_control()), collapse = ', '), "\n")
      }
    } else if (!use.grad && tree.unbiased) { 
      if (!setequal(names(mob_control()), names(tree.control))) {
        stop("Argument 'tree.control' should be a list containing named elements ", paste(names(mob_control()), collapse = ', '), "\n")
      }
    } else if (!tree.unbiased) {
      if (!setequal(names(rpart.control()), names(tree.control))) {
        stop("Argument 'tree.control' should be a list containing names elements ", paste(names(rpart.control()), collapse = ', '), "\n")
      }
    }
    if (use.grad) { ## if ctree or rpart are employed:
      tree.control$maxdepth <- maxdepth[1L]      
    } else if (tree.unbiased) { ## if glmtree is employed:
      tree.control$maxdepth <- maxdepth[1L] + 1L
    }
    if (tree.unbiased) {
      tree.control$mtry <- mtry
    } else if (mtry != Inf) {
      warning("Argument 'tree.unbiased' was set to FALSE, so rpart is employed for tree induction, and value specified for 'mtry' will be ignored.\n")
      mtry <- Inf
    }
  }
  
  if (!tree.unbiased && !use.grad && learnrate > 0) {
    stop("Employing the rpart algorithm with a learnrate > 0 without gradient boosting is not supported.\n")
  }
  
  
  ######################################
  ## Prepare data, formula and family ##
  ######################################

  ## prepare model frame:
  data <- model.frame(Formula::as.Formula(formula), data = data, na.action = NULL)

  ## Coerce character and logical variables to factors:
  if (any(char_names <- sapply(data, is.character))) {
    char_names <- names(data)[char_names]
    warning("The following variables were of class 'character' and will be coerced to 'factor': ", paste(char_names, collapse = " "), "\n")
    data[ , char_names] <- sapply(data[ , char_names], factor)
  }
  if (any(logic_names <- sapply(data, is.logical))) {
    logic_names <- names(data)[logic_names]
    warning("The following variables were of class 'logical' and will be coerced to 'factor': ", paste(logic_names, collapse = " "), "\n")
    data[ , logic_names] <- sapply(data[ , logic_names], factor)
  } 
  
  ## Coerce ordered categorical variables to numeric:
  if (ordinal) {
    if (any(ord_var_inds <- sapply(data, is.ordered))) {
      data[ , ord_var_inds] <- sapply(data[ , ord_var_inds], as.numeric)
    }
  }
  
  ## get response variable name(s):
  y_names <- names(data)[attr(attr(data, "terms"), "response")]
  if (family == "mgaussian" || length(y_names) == 0) {
    y_names <- attr(terms(Formula(formula), rhs = 0, data = data), "term.labels")
    family <- "mgaussian"
  }
  
  ## get predictor variable names:
  if (family == "cox" || is.Surv(data[, y_names])) {
    x_names <- attr(attr(data, "terms"), "term.labels")
  } else if (use_glmertree) {
    ## TODO: This does, but should not, remove all functions used in formula:
    x_names <- all.vars(formula[[3L]][[3L]])
  } else {
    x_names <- attr(terms(Formula(formula), lhs = 0, data = data), "term.labels")
  }
  
  ## expand dot and put ticks around variables within functions, if present:
  if (!(use_glmertree || family == "mgaussian")) {
    formula <- formula(data)
  } else if (family == "mgaussian") {
    formula <- Formula(formula(paste0(
      paste0(paste0("`", y_names, "`"), collapse = " + "), 
      " ~ ", 
      paste0(paste0("`", x_names, "`"), collapse = " + "))))
  }

  ## get sample size:
  n <- nrow(data)

  ## check and set correct family:
  if (is.null(cl$family)) {
    if (length(y_names) == 1L) {
      if (is.factor(data[,y_names])) { # then family should be bi- or multinomial
        if (is.ordered(data[,y_names])) {
          warning("An ordered factor was specified as the response variable, which will be treated as an unordered factor response.")
          data[,y_names] <- factor(data[,y_names], ordered = FALSE)
        } 
        if (nlevels(data[,y_names]) == 2L) {
            family <- "binomial"
        } else if (nlevels(data[,y_names]) > 2L) {
            family <- "multinomial"
          }
      } else if (is.Surv(data[,y_names])) { # then family should be cox
        family <- "cox"
      } else if (!is.numeric(data[,y_names])) { # then response is not a factor, survival or numeric
        warning("The response variable specified through argument 'formula' should be of class numeric, factor or Surv.")
      }
    } else if (length(y_names) > 1L) { # multiple responses specified, should be numeric
      if (all(sapply(data[,y_names], is.numeric))) {
        family <- "mgaussian"
      } else {
        warning("Multiple response variables were specified, but not all were (but should be) numeric.\n")
      }
    }
    
  } else { # family was specified, check if correct;
    
    if (family[1L] == "gaussian") {
      if (length(y_names) > 1L) {
        warning("Argument 'family' was set to 'gaussian', but multiple response variables were specified in 'formula'.\n")        
      }
      if (!is.numeric(data[,y_names])) { # then family should be poisson or gaussian
        warning("Argument 'family' was set to 'gaussian', but the response variable specified in 'formula' is not of class numeric.\n")
      }
    } else if (family[1L] == "poisson") {
      if (length(y_names) > 1L) {
        warning("Argument 'family' was set to 'poisson', but multiple response variables were specified, which is not supported.\n")
      }
      if (!isTRUE(all.equal(round(data[,y_names]), data[,y_names]))) {
        warning("Argument 'family' was set to 'poisson', but the response variable specified in 'formula' is non-integer.\n")
      }
    } else if (family[1L] == "binomial") {
      if (length(y_names) > 1L) {
        warning("Argument 'family' was set to 'binomial', but multiple response variables were specified, which is not supported.\n")
      } else if (!is.factor(data[,y_names])) {
        warning("Argument 'family' was set to 'binomial', but the response variable specified is not a factor.\n")
      } else if (nlevels(data[,y_names]) != 2L) {
        warning("Argument 'family' was set to 'binomial', but the response variable has ", nlevels(data[,y_names]), " levels.\n")        
      }
    } else if (family[1L] == "multinomial") {
      if (length(y_names) > 1L) {
        warning("Argument 'family' was set to 'multinomial', but multiple response variables were specified, which is not supported.\n")
      } else if (!is.factor(data[,y_names])) {
        warning("Argument 'family' was set to 'multinomial', but the response variable specified is not a factor.\n")
      } else if (nlevels(data[,y_names]) < 3L) {
        warning("Argument 'family' was set to 'multinomial', but the response variable has ", nlevels(data[,y_names]), " levels.\n")
      }
    } else if (family[1L] == "cox") {
      if (length(y_names) > 1L) {
        warning("Argument 'family' was set to 'cox', but multiple response variables were specified, which is not supported.\n")
      } else if (!is.Surv(data[ , y_names])) {
        warning("Argument 'family' was set to 'cox', but the response variable specified is not of class Surv.\n")
      }
    } else if (family == "mgaussian") {
      if (length(y_names) == 1L) {
        warning("Argument 'family' was set to 'mgaussian', but only a single response variable was specified.\n")
      } else if (!all(sapply(data[,y_names], is.numeric))) {
        warning("Argument 'family' was set to 'mgaussian', but not all response variables specified are numeric.\n")
      }
    }
  }



  ## Check specification of tree growing algorithms employed:
  if (!tree.unbiased) { # rpart is employed
    if (family == "mgaussian") {
      stop("Employing rpart algorithm for rule induction with a multivariate response variable is not supported. Set argument 'tree.unbiased' to TRUE and argument 'use.grad' to FALSE.\n")
    } else if (learnrate > 0 && family == "multinomial") {
      stop("Employing rpart algorithm for rule induction with a multinomial response variable and learnrate > 0 is not supported. Set argument 'learnrate' to 0, or arguments 'tree.unbiased' and 'use.grad' to TRUE.\n")
    }
  } else if (!use.grad) { # (g)lmtree is employed
    if (family == "multinomial") {
      stop("Employing (g)lmtree for rule induction with a multinomial response variable is not supported. Set argument 'use.grad' to TRUE for multivariate responses.\n")
    } else if (family == "mgaussian") {
      stop("Employing (g)lmtree for rule induction with a multivariate response variable is not supported. Set argument 'use.grad' to TRUE for multivariate responses.\n")
    } else if (family == "cox") {
      stop("Employing (g)lmertree for rule induction with a survival response is not supported. Set argument 'use.grad' to TRUE for a survival response.\n")
    }
  }
  
  if (family == "cox") {
    if (!requireNamespace("survival", quietly = TRUE)) {
      stop("For fitting a prediction rule ensemble with a survival response, package survival should be installed and loaded.\n")    
    }
    if (learnrate > 0) {
      if (!requireNamespace("mboost", quietly = TRUE)) {
        stop("For fitting a prediction rule ensemble with a survival response and learning rate > 0, package mboost should be installed.\n")
      }
    }
  }

  ## Prevent response from being interpreted as count by ctree or rpart:
  if (learnrate == 0 && family == "gaussian" && (!(tree.unbiased && !use.grad))) { # if glmtree is not employed
    if (isTRUE(all.equal(round(data[, y_names]), data[, y_names]))) { # if response passes integer test
      data[, y_names] <- data[, y_names] + 0.01 # add small constant to response to prevent response being interpreted as count by ctree or rpart
      small_constant_added <- 0.01
    } else {
      small_constant_added <- FALSE
    }
  } else {
    small_constant_added <- FALSE
  }

  if (any(is.na(data))) {
    weights <- weights[complete.cases(data)]
    data <- data[complete.cases(data),]
    n <- nrow(data)
    warning("Some observations have missing values and have been removed from the data. New sample size is ", n, ".\n", immediate. = TRUE)
  }

  if (verbose) {
    if (family == "gaussian") {
      cat("A rule ensemble for prediction of a continuous response will be created.\n")
    } else if (family == "poisson") {     
      cat("A rule ensemble for prediction of a count response will be created.\n")
    } else if (family == "binomial") {
      cat("A rule ensemble for prediction of a binary categorical response will be created.\n")
    } else if (family == "multinomial") {
      cat("A rule ensemble for prediction of a categorical response with > 2 levels will be created.\n")
    } else if (family == "mgaussian") {
      cat("A rule ensemble for prediction of a multivariate continuous response will be created.\n")
    } else if (family == "cox") {
      cat("A rule ensemble for prediction of a survival response will be created.\n")
    }
  }
  
    
  #############################
  ## Derive prediction rules ##
  #############################
  
  if (type == "linear") {
    rules <- NULL
    rulevars <- NULL
    
  } else {
    if (use_glmertree) {
      rule_object <- pre_rules_mixed_effects(formula = formula, 
                                data = data,
                                y_names = y_names,
                                x_names = x_names,
                                learnrate = learnrate, 
                                par.init = par.init, 
                                sampfrac = sampfrac, 
                                mtry = mtry,
                                weights = weights, 
                                ntrees = ntrees, 
                                tree.control = tree.control, 
                                maxdepth = maxdepth,
                                use.grad = use.grad, 
                                family = family, 
                                verbose = verbose, 
                                removeduplicates = removeduplicates, 
                                removecomplements = removecomplements)
    } else {
      rule_object <- pre_rules(formula = formula, 
                               data = data,
                               weights = weights,
                               y_names = y_names,
                               x_names = x_names,
                               learnrate = learnrate, 
                               par.init = par.init, 
                               sampfrac = sampfrac, 
                               mtry = mtry,
                               maxdepth = maxdepth,
                               ntrees = ntrees, 
                               tree.control = tree.control, 
                               use.grad = use.grad, 
                               family = family, 
                               verbose = verbose, 
                               removeduplicates = removeduplicates, 
                               removecomplements = removecomplements,
                               tree.unbiased = tree.unbiased,
                               return.dupl.compl = TRUE, 
                               sparse = sparse)
    }
    rules <- rule_object$rules
    rulevars <- rule_object$rulevars
    
  }

  #########################################################################
  ## Prepare rules, linear terms, outcome variable and perform selection ##
  #########################################################################
  
  if (type == "rules" && length(rules) == 0) {
    warning("No prediction rules could be derived from the data.")
    return(NULL)
  }
  
  if (is.numeric(small_constant_added)) {
    data[, y_names] <- data[, y_names] - small_constant_added
  }
  
  ## Prepare right formula if glmertree was used for tree induction
  if (use_glmertree) {
    modmat_f <- as.formula(paste(formula[[2]], formula[[1]], paste(x_names, collapse = "+")))
  } else {
    modmat_f <- formula
  }
  modmat_data <- get_modmat(
    formula = modmat_f, 
    data = data, 
    rules = rules, 
    type = type, 
    winsfrac = winsfrac, 
    x_names = x_names,
    y_names = y_names,
    normalize = normalize, 
    sparse = sparse, 
    rulevars = rulevars,
    confirmatory = confirmatory)
  rm(rulevars)
  
  x_scales <- modmat_data$x_scales
  wins_points <- modmat_data$wins_points
  
  
  ############################
  ### Fit regression model ###
  ############################
  
  ### TODO: allow for forward selection:
  ## Include additional argument regression = "glmnet"
  ## which also takes argument "stepAIC"
  ## then number of terms (i.e., steps argument) should be specified
  ## But would be nice to always take e.g., 100 steps, 
  ## and then select number of terms with penalty.par.val in print etc. 
  ## would allow only for ""continuous
  
  #if (regression == "stepAIC") {
  #  if (family %in% c("gaussian", "binomial")) {
  #    data <- cbind(modmat_data$y, modmat_data$x)
  #    lm_full <- lm()
  #    lm_intercept <- lm()
  #    MASS::stepAIC(object = lm_intercept, scope = list(upper = lm_full, lower = lm_intercept),
  #                  direction = "forward", type)    
  #    ## with stepAIC seems tricky to save intermediate models.
  #    ## Create a loop with MASS::addterm
  #    
  #    
  #  }
  #  
  #
  #} else if (regression == "glmnet") {
  
  y <- modmat_data$y
  x <- modmat_data$x  

  if (is.null(cl$penalty.factor)) {
    penalty.factor <- rep(1L, times = ncol(x))
  } 
  if (!is.null(confirmatory)) {
    penalty.factor[which(colnames(x) %in% confirmatory)] <- 0L
  }
  
  # check whether there's duplicates in the variable names:
  # (can happen, for example, due to labeling of dummy indicators for factors)
  if (!(length(unique(colnames(x))) == length(colnames(x)))) { 
    warning("There are variables in the model with overlapping variable names. If predictor variables of type factor were specified with numbers in their name, consider renaming these and rerunning the analysis. See 'Details' under ?pre.") 
  } 
  glmnet.fit <- cv.glmnet(x, y, nfolds = nfolds, weights = weights, 
                          family = family, parallel = par.final, 
                          standardize = standardize, 
                          penalty.factor = penalty.factor,
                          ...)
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
  #}
  
  
  ####################
  ## Return results ##
  ####################

  result <- list(glmnet.fit = glmnet.fit, call = cl, weights = weights, 
                 data = data, normalize = normalize, x_scales = x_scales, 
                 type = type, x_names = x_names, y_names = y_names, 
                 modmat = x, wins_points = wins_points, family = family, 
                 formula = formula)
  if (type != "linear" & length(rules) > 0) {
    result$complements.removed <- rule_object$complements.removed
    result$duplicates.removed <- rule_object$duplicates.removed
    result$rules <- data.frame(rule = names(rules), 
                               description = rules, 
                               stringsAsFactors = FALSE)
  } 
  
  class(result) <- "pre"
  return(result)
}



get_modmat <- function(
  # Pass these if you already have an object
  wins_points = NULL, x_scales = NULL,
  # These should be passed in all calls
  formula, data, rules, type, x_names, winsfrac, normalize, 
  # Response variable names is optional:
  y_names = NULL, sparse = FALSE, rulevars = NULL,
  confirmatory = NULL) {
  
  ## TODO: is next line necessary? Only for converting ordered categorical vars to linear terms?
  ## Is only used twice, in get_rulemat function below.
  ## Perhaps should be conditional on ordinal argument of pre()?
  data_org <- data # needed to evaluate rules later
  
  # convert ordered categorical predictor variables to linear terms:
  data[,sapply(data, is.ordered)] <- 
    as.numeric(data[,sapply(data, is.ordered)])
  
  ## Add confirmatory terms:
  if (!is.null(confirmatory)) {
    ## Get confirmatory rules and variables:
    conf_vars <- confirmatory[confirmatory %in% x_names]
    conf_rules <- confirmatory[!(confirmatory %in% conf_vars)]
    if (length(conf_rules) > 0L) {
      ## Evaluate rules (add to model matrix later):
      eval_rules <- function(data, rule_descr) {
        1L * with(data, eval(parse(text = rule_descr)))
      }
      conf_rules <- mapply(FUN = eval_rules, rule_descr = conf_rules, 
                           MoreArgs = list(data = data))
    }
  }
  
  #####
  # Perform winsorizing and normalizing
  if (type != "rules" && any(sapply(data[,x_names], is.numeric))) {
    #####
    # if type is not rules, linear terms should be prepared:
    
    # Winsorize numeric variables (section 5 of F&P(2008)):
    if (winsfrac > 0) {
      
      if (miss_wins_points <- is.null(wins_points)) {
        wins_points <- data.frame(varname = x_names, value = NA, lb = NA, ub = NA)
      }
      
      j <- 0L
      tol <- sqrt(.Machine$double.eps)
      for (i in x_names) {
        j <- j + 1L
        if (is.numeric(data[[i]])) {
          if (miss_wins_points) {
            lim <- quantile(data[, i], probs = c(winsfrac, 1 - winsfrac))
            wins_points$value[j] <- paste(lim[1L], "<=", i, "<=", lim[2L])
            lb <- lim[1L]
            ub <- lim[2L]
            if (ub - lb < tol) {
              ## If lower and upper bound are equal, do not winsorize and issue warning:
              warning("Variable ", x_names[j], " will be winsozired employing winsfrac = 0, to prevent reducing the variance of its linear term to 0.", immediate. = TRUE)
              wins_points$lb[j] <- min(data[ , i])
              wins_points$ub[j] <- max(data[ , i])
            } else {
              wins_points$lb[j] <- lb
              wins_points$ub[j] <- ub
              data[ , i][data[ , i] < lb] <- lb
              data[ , i][data[ , i] > ub] <- ub
            }
          } else {
            data[ , i][data[ , i] < wins_points$lb[j]] <- wins_points$lb[j]
            data[ , i][data[ , i] > wins_points$ub[j]] <- wins_points$ub[j]
          }
        }
      }
    }
    
    # normalize numeric variables:
    if (normalize) { 
      # Normalize linear terms (section 5 of F&P08), if there are any:
      needs_scaling <- x_names[sapply(data[ , x_names], is.numeric)]
      if (length(needs_scaling) > 0) {
        if (is.null(x_scales)) {
          x_scales <- apply(
            data[, needs_scaling, drop = FALSE], 2L, sd, na.rm = TRUE) / 0.4
        }
        ## check if variables have zero variance (if so, do not scale):
        tol <- sqrt(.Machine$double.eps)
        almost_zero_var_inds <- which(x_scales < tol)
        if (length(almost_zero_var_inds) > 0) {
          # print warning and set all those x_scales to 1
          warning("Variable(s) ", needs_scaling[almost_zero_var_inds], " have sd < ", tol, " and will not be normalized.")  
          # omit from needs_scaling:
          x_scales[almost_zero_var_inds] <- 1
        }
        data[ , needs_scaling] <- scale(
          data[ , needs_scaling, drop = FALSE], center = FALSE, scale = x_scales)
      }
    }
  }

  ## Combine rules and variables:
  if (length(rules) == 0L) rules <- NULL
  if (type == "linear" || is.null(rules)) {
    if (sparse) {
      mf <- model.frame(Formula(formula), data)
      x <- model.Matrix(terms(mf), mf, sparse = TRUE)
      rm(mf)
      
    } else {
      x <- model.matrix(Formula(formula), data = data)
    }
  } else if (type %in% c("both", "rules") && !is.null(rules)) {
    if (sparse) {
      x <- if (is.null(rulevars)) { 
        .get_rules_mat_sparse(data_org, rules) 
        } else {
          rulevars
        }
      if (type == "both") {
        mf <- model.frame(Formula(formula), data)
        x <- cbind(model.Matrix(terms(mf), mf, sparse = TRUE), x)
        rm(mf)
      }
    } else { 
      x <- if (is.null(rulevars)) {
        .get_rules_mat_dense(data_org, rules)
        } else { 
          rulevars
        }
      if (type == "both") {
        x <- cbind(model.matrix(Formula(formula), data = data), x)
      }
    }
    
  } else { 
    stop("not implemented with type ", sQuote(type), " and is.null(rules) is ", 
         sQuote(is.null(rules)))
  }
  
  #####
  # Remove intercept
  x <- x[, colnames(x) != "(Intercept)", drop = FALSE]

  ## Add confirmatory terms:
  if (!is.null(confirmatory)) {
    ## Get confirmatory rules and variables:
    #conf_vars <- confirmatory[confirmatory %in% x_names]
    #conf_rules <- confirmatory[!(confirmatory %in% conf_vars)]
    if (length(conf_vars) > 0L) {
      ## If type == "rules", add variables to model matrix:
      if (type == "rules") {
        x <- cbind(x, as.matrix(data[conf_vars], rownames.force = FALSE))
        ## TODO: Then check if it should be winsorized or normalized, or issue warning if not
      }
    }
    if (is.matrix(conf_rules)) {
      if (!sparse) {
        x <- cbind(x, conf_rules)
      } else {
        ## TODO: implement something for sparse matrices.
      }
    }
  }
  
  if (is.null(y_names)) {
    y <- NULL
  } else {
    y <- data[ , y_names]
    if (is.Surv(y) || length(y_names) > 1L) {
      y <- as.matrix(y)
    }
  }
  
  if (!exists("wins_points", inherits = FALSE)) { wins_points <- NULL }
  if (!exists("x_scales", inherits = FALSE)) { x_scales <- NULL }
  
  list(x = x, y = y, x_scales = x_scales, wins_points = wins_points)
}

.get_rules_mat_dense <- function(data, rules){
  if(length(rules) == 0)
    return(NULL)
  
  expr <- parse(text = paste0("cbind(", paste0(rules, collapse = ", "), ")"))
  x <- eval(expr, data)
  colnames(x) <- names(rules)
  x
}

.get_rules_mat_sparse <- function(data, rules){
  if(length(rules) == 0)
    return(NULL)
  
  # See https://stackoverflow.com/a/8844057/5861244. 
  #
  # if all rules where binary then we could use the `lsparseMatrix-classes`. 
  # However, this will not work when we call `glmnet` as it requires a double
  # matrix`
  expr <- paste0("cbind_sparse_vec(", paste0(
    'as(as.numeric(', rules, '), "sparseVector")', collapse = ", "), ")")
  x <- eval(parse(text = expr), data)
  colnames(x) <- names(rules)
  x
}




## Rule learner for pre:
pre_rules <- function(formula, data, weights = rep(1, nrow(data)),
                      y_names, x_names, 
                      learnrate = .01, par.init = FALSE, sampfrac = .5, 
                      mtry = Inf, maxdepth = 3L, ntrees = 500, 
                      tree.control = ctree_control(), use.grad = TRUE, 
                      family = "gaussian", verbose = FALSE, 
                      removeduplicates = TRUE, removecomplements = TRUE,
                      tree.unbiased = TRUE, return.dupl.compl = FALSE, 
                      sparse = FALSE) {
  
  n <- nrow(data)
  
  ## Prepare glmtree arguments, if necessary:
  if (!use.grad && tree.unbiased) {
    glmtree_args <- mob_control(maxdepth = maxdepth[1] + 1, mtry = mtry)
    glmtree_args$formula <- formula(paste(paste(y_names, " ~ 1 |"), 
                                          paste(x_names, collapse = "+")))
    if (!family == "gaussian") {
      glmtree_args$family <- family      
    }
  } else {
    glmtree_args <- NULL
  }

  ## Set up subsamples (outside of the loop):
  subsample <- list()
  for (i in 1:ntrees) {
    if (is.function(sampfrac)) {
      subsample[[i]] <- sampfrac(n = n, weights = weights)
    }
    else if (sampfrac == 1) {
      subsample[[i]] <- sample(1:n, size = n, replace = TRUE, 
                               prob = weights)
    }
    else if (sampfrac < 1) {
      subsample[[i]] <- sample(1:n, size = round(sampfrac * n), 
                               replace = FALSE, prob = weights)
    }
  }  
  
  ## Grow trees:
  if (learnrate == 0) {
    
    ## Set up rule learning function:
    fit_tree_return_rules <- function(formula, data, family = NULL, 
                                      use.grad = TRUE, tree.unbiased = TRUE,
                                      glmtree_args = NULL, tree.control = NULL) {
      if (tree.unbiased) {
        if (use.grad) { # employ ctree
          tree <- ctree(formula = formula, data = data, control = tree.control)
          return(list.rules(tree, removecomplements = removecomplements))
        } else { # employ (g)lmtree
          glmtree_args$data <- data
          if (family == "gaussian") {
            tree <- do.call(lmtree, args = glmtree_args)
          } else {
            tree <- do.call(glmtree, args = glmtree_args)
          }
          return(list.rules(tree, removecomplements = removecomplements))
        }
      } else { # employ rpart
        tree <- rpart(formula = formula, data = data, control = tree.control)
        paths <- path.rpart(tree, nodes = rownames(tree$frame), print.it = FALSE)
        paths <- unname(sapply(sapply(paths, `[`, index = -1), paste, collapse = " & ")[-1])
        if (removecomplements) {
          ## Omit first rule, as it is the complement of a later rule, by definition:
          paths <- paths[-1]
        }
        return(paths)
      }
    }
    
    if (par.init) { # compute in parallel:
      rules <- foreach::foreach(i = 1:ntrees, .combine = "c", .packages = c("partykit", "pre")) %dopar% {
        
        if (length(maxdepth) > 1L) {
          if (use.grad) {
            tree.control$maxdepth <- maxdepth[i]
          } else if (tree.unbiased) {
            glmtree_args$maxdepth <- maxdepth[i] + 1L
          }
        }
        fit_tree_return_rules(data[subsample[[i]], ], 
                              formula = formula, 
                              family = family, 
                              use.grad = use.grad, 
                              tree.unbiased = tree.unbiased, 
                              glmtree_args = glmtree_args, 
                              tree.control = tree.control)
      }
      
    } else { # compute in serial:
      
      rules <- c()
      for (i in 1:ntrees) {
        
        if (length(maxdepth) > 1L) {
          if (use.grad) {
            tree.control$maxdepth <- maxdepth[i]
          } else if (tree.unbiased) {
            glmtree_args$maxdepth <- maxdepth[i] + 1L
          }
        }
        rules <- c(rules, 
                   fit_tree_return_rules(data[subsample[[i]], ], 
                                         formula = formula, 
                                         family = family, 
                                         use.grad = use.grad, 
                                         tree.unbiased = tree.unbiased, 
                                         glmtree_args = glmtree_args, 
                                         tree.control = tree.control))
        
      }
    }
  
  } else { # learnrate > 0:
    
    rules <- c() # initialize with empty rule vector
    
    if (use.grad) { ## use ctrees or rpart with y_learn and eta:
      
      data_with_y_learn <- data
      ## set initial y and eta value:
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
        ## create dummy variables:
        y <- model.matrix(as.formula(paste0(" ~ ", y_names, " - 1")), data = y)
        ## adjust formula used by ctree to involve multiple response variables:
        formula <- as.formula(paste(paste(colnames(y), collapse = " + "), "~", 
                                    paste(x_names, collapse = " + ")))
        ## get y_learn:
        eta_0 <- get_intercept_multinomial(y, weights)
        eta <- t(replicate(n = nrow(y), expr = eta_0))
        p_0 <- 1 / (1 + exp(-eta))
        for (i in 1:ncol(y)) {
          y[,i] <- ifelse(y[,i] == 1, log(p_0[,i]), log(1 - p_0[,i]))
        }
        ## omit original response and include dummy-coded response in data:
        data_with_y_learn <- cbind(data[,-which(names(data)== y_names)], y)
        multinomial_y_names <- names(y)
      } else if (family == "mgaussian") {
        y <- data[,y_names]
        eta_0 <- apply(y, 2, weighted.mean, weights = rep(1, nrow(y)))
        eta <- t(replicate(n = nrow(y), expr = eta_0))
        data_with_y_learn[,y_names] <- y - eta
      } else if (family == "cox") {
        ## Adjust formula used by ctree and rpart:
        formula <- as.formula(paste0("pseudo_y ~ ", 
                                     paste0(x_names, collapse = " + ")))
        y <- data[,y_names]
        eta_0 <- 0
        eta <- rep(0, times = nrow(data))
        ngradient_CoxPH <- mboost::CoxPH()@ngradient
        ## omit original response and include pseudo-y:
        data_with_y_learn <- cbind(data[,-which(names(data)== y_names)], y)
        data_with_y_learn$pseudo_y <- ngradient_CoxPH(y = y, f = eta, w = weights)
      }

      for(i in 1:ntrees) {

        if (length(maxdepth) > 1L) {
          tree.control$maxdepth <- maxdepth[i]
        }
        # Grow tree on subsample:
        if (tree.unbiased) {
          tree <- ctree(formula, control = tree.control,
                        data = data_with_y_learn[subsample[[i]], ])
          # Collect rules:
          rules <- c(rules, list.rules(tree, removecomplements = removecomplements))
        } else {
          tree <- rpart(formula, control = tree.control,
                        data = data_with_y_learn[subsample[[i]], ])
          paths <- path.rpart(tree, nodes = rownames(tree$frame), print.it = FALSE, pretty = 0)
          paths <- unname(sapply(sapply(paths, `[`, index = -1), paste, collapse = " & ")[-1])
          if (removecomplements) {
            paths <- paths[-1]
          }
          
          rules <- c(rules, paths)
        }
        
        ## Update eta and y_learn:
        eta <- eta + learnrate * predict(tree, newdata = data)
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
      
    } else { ## use.grad is FALSE, employ (g)lmtrees with offset:
      
      ## initialize with 0 offset:
      offset <- rep(0, times = nrow(data))
      
      for(i in 1:ntrees) {
        
        # Take subsample of dataset:
        glmtree_args$data <- data[subsample[[i]],]
        glmtree_args$offset <- offset[subsample[[i]]] 
        if (length(maxdepth) > 1L) {
          glmtree_args$maxdepth <- maxdepth[i] + 1L
        }
        # Grow tree on subsample:
        if (family == "gaussian") {
          tree <- do.call(lmtree, args = glmtree_args)      
        } else {
          tree <- do.call(glmtree, args = glmtree_args) 
        }
        # Collect rules:
        rules <- c(rules, list.rules(tree, removecomplements = removecomplements))
        # Update offset (note: do not use a dataset which includes the offset for prediction!!!):
        if (learnrate > 0) {
          if (family == "gaussian") {
            offset <- offset + learnrate * 
              suppressWarnings(predict(tree, newdata = data))
          } else {
            offset <- offset + learnrate * 
              suppressWarnings(predict(tree, newdata = data, type = "link"))
          }
        }
      }
    }
  }
  
  if (length(rules) > 0) {
    # Keep unique, non-empty rules only:
    rules <- unique(rules[!rules==""])
    if (sparse) {
      rules <- .get_most_sparse_rule(rules, data)
    }
    ## Adjust rule format of rpart rules:
    if (!tree.unbiased) {
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
  
    #complements.removed <- rules_obj$complements.removed
    #duplicates.removed <- rules_obj$duplicates.removed
    #rules <- rules_obj$rules
    #rulevars <- rules_obj$rulevars
  
    if (verbose && (removeduplicates || removecomplements)) 
      cat("\n\nA total of", length(rules_obj$duplicates.removed) + length(rules_obj$complements.removed), "generated rules were perfectly collinear with earlier rules and removed from the initial ensemble. \n($duplicates.removed and $complements.removed show which, if any).")
    
    if (verbose)
      cat("\n\nAn initial ensemble consisting of", length(rules_obj$rules), "rules was successfully created.")  
    
  } else {
    warning("No prediction rules could be derived from dataset.", immediate. = TRUE)  
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







pre_rules_mixed_effects <- function(formula, data, family = "gaussian", 
                                    y_names, x_names, learnrate = .01, 
                                    sampfrac = .5, 
                                    weights = rep(1, nrow(data)), 
                                    mtry = Inf, maxdepth = 3L, ntrees = 500,
                                    tree.control = ctree_control(mtry = mtry, maxdepth = maxdepth[1]), 
                                    use.grad = TRUE, verbose = FALSE, 
                                    removeduplicates = TRUE, 
                                    removecomplements = TRUE, 
                                    par.init = FALSE, sparse = FALSE) {
  
  n <- nrow(data)
  
  ## Prepare arguments:
  glmertree_args <- tree.control    
  glmertree_args$formula <- formula
  if (family != "gaussian") {
    glmertree_args$family <- family      
  }
  
  ## TODO: Allow for boosting, maybe allow for applying learnrate to the mixed- or only fixed-effects predictions
  
  ## Setup samples:
  subsample <- list()
  for (i in 1:ntrees) {
    if(is.function(sampfrac)) {
      subsample[[i]] <- sampfrac(n = n, weights = weights) 
    } else if (sampfrac == 1) { # then bootstrap
      subsample[[i]] <- sample(1:n, size = n, replace = TRUE, prob = weights)
    } else if (sampfrac < 1) { # else subsample
      subsample[[i]] <- sample(1:n, size = sampfrac*n, replace = FALSE, prob = weights)
    }
  }
  
  rules <- c()
  
  if (par.init) { # compute in parallel:
    
    rules <- foreach::foreach(i = 1:ntrees, .combine = "c", .packages = c("partykit", "pre")) %dopar% {
      
      # Prepare call:
      glmertree_args$data <- data[subsample[[i]], ]
      if (length(maxdepth) > 1) {
        glmertree_args$maxdepth <- maxdepth[i] + 1        
      }
      # Fit tree:
      if (family == "gaussian") {
        tree <- do.call(glmertree::lmertree, args = glmertree_args)$tree
      } else {
        tree <- do.call(glmertree::glmertree, args = glmertree_args)$tree
      }
      # Collect rules:
      list.rules(tree, removecomplements = removecomplements)
        
    } 
  } else { # do not compute in parallel:
      
    rules <- c()
    for (i in 1:ntrees) {

      # prepare call:
      glmertree_args$data <- data[subsample[[i]], ]
      if (length(maxdepth) > 1) {
        glmertree_args$maxdepth <- maxdepth[i] + 1        
      }
      # Fit tree:
      if (family == "gaussian") {
        tree <- do.call(glmertree::lmertree, args = glmertree_args)$tree
      } else {
        tree <- do.call(glmertree::glmertree, args = glmertree_args)$tree
      }
      # Collect rules:
      rules <- c(rules, list.rules(tree, removecomplements = removecomplements))
    
    }
  } 
    
  # Keep unique, non-empty rules only:
  rules <- unique(rules[!rules==""])
  if (verbose) {
    cat("\nA total of", ntrees, "trees and ", length(rules), "rules were generated initially.")
  }
  
  rules_obj <- delete_duplicates_complements(
    rules = rules, data = data, 
    removecomplements = removecomplements, 
    removeduplicates = removeduplicates, 
    return.dupl.compl = TRUE, sparse = sparse, keep_rulevars = TRUE)
  
  complements.removed <- rules_obj$complements.removed
  duplicates.removed <- rules_obj$duplicates.removed
  rules <- rules_obj$rules
  rulevars <- rules_obj$rulevars
  
  if (verbose && (removeduplicates || removecomplements)) 
    cat("\n\nA total of", length(duplicates.removed) + length(complements.removed), "generated rules were perfectly collinear with earlier rules and removed from the initial ensemble. \n($duplicates.removed and $complements.removed show which, if any).")
  
  if (verbose) 
    cat("\n\nAn initial ensemble consisting of", length(rules), "rules was succesfully created.")  

  # check if any rules were generated:
  if (length(rules) == 0)
    warning("No prediction rules could be derived from dataset.", immediate. = TRUE)
  
  rules_obj
}


#' Sampling function generator for specifyinf varying maximum tree depth 
#' in a prediction rule ensemble (pre)
#' 
#' \code{maxdepth_sampler} generates a random sampling function, governed
#' by a pre-specified average tree depth.
#' 
#' @param av.no.term.nodes integer of length one. Specifies the average 
#' number of terminal nodes in trees used for rule inducation.
#' @param av.tree.depth integer of length one. Specifies the average maximum
#' tree depth in trees used for rule induction.
#' @return Returns a random sampling function with single argument 'ntrees',
#' which can be supplied to the \code{maxdepth} argument of function 
#' \code{\link{pre}} to specify varying tree depths.
#' @details The original RuleFit implementation varying tree sizes for
#' rule induction. Furthermore, it defined tree size in terms of the number
#' of terminal nodes. In contrast, function \code{\link{pre}} defines the 
#' maximum tree size in terms of a (constant) tree depth. Function 
#' \code{maxdepth_sampler} allows for mimicing the behavior of the
#' orignal RuleFit implementation. In effect, the maximum tree depth is 
#' sampled from an exponential distribution with learning rate 
#' \eqn{1/(\bar{L}-2)}, where \eqn{\bar{L} \ge 2} represents the
#' average number of terminal nodes for trees in the ensemble. See
#' Friedman & Popescu (2008, section 3.3). 
#' @references Friedman, J. H., & Popescu, B. E. (2008). Predictive learning 
#' via rule ensembles. \emph{The Annals of Applied Statistics, 2}(3), 916-954.
#' @seealso \code{\link{pre}}
#' @examples
#' ## RuleFit default is max. 4 terminal nodes, on average:
#' func1 <- maxdepth_sampler()
#' set.seed(42)
#' func1(10)
#' mean(func1(1000))
#' 
#' ## Max. 16 terminal nodes, on average (equals average maxdepth of 4):
#' func2 <- maxdepth_sampler(av.no.term.nodes = 16L)
#' set.seed(42)
#' func2(10)
#' mean(func2(1000))
#' 
#' ## Max. tree depth of 3, on average:
#' func3 <- maxdepth_sampler(av.tree.depth = 3)
#' set.seed(42)
#' func3(10)
#' mean(func3(1000))
#' 
#' ## Max. 2 of terminal nodes, on average (always yields maxdepth of 1):
#' func4 <- maxdepth_sampler(av.no.term.nodes = 2L)
#' set.seed(42)
#' func4(10)
#' mean(func4(1000))
#' 
#' \donttest{## Create rule ensemble with varying maxdepth:
#' set.seed(42)
#' airq.ens <- pre(Ozone ~ ., data = airquality[complete.cases(airquality),],
#'                 maxdepth = func1)
#' airq.ens}
#' @export
maxdepth_sampler <- function(av.no.term.nodes = 4L, av.tree.depth = NULL) {
  function(ntrees, ...) {
    if (!is.null(av.tree.depth)) {
      av.no.term.nodes <- 2^av.tree.depth
    }
    ceiling(log(2 + floor(rexp(ntrees, rate = 1 / (av.no.term.nodes - 2))), base = 2))
  }
}



#' Print method for objects of class pre
#'
#' \code{print.pre} prints information about the generated prediction rule 
#' ensemble to the command line
#' 
#' @param x An object of class \code{\link{pre}}.
#' @param penalty.par.val character or numeric. Value of the penalty parameter
#' \eqn{\lambda} to be employed for selecting the final ensemble. The default
#' \code{"lambda.min"} employs the \eqn{\lambda} value within 1 standard
#' error of the minimum cross-validated error. Alternatively, 
#' \code{"lambda.min"} may be specified, to employ the \eqn{\lambda} value
#' with minimum cross-validated error, or a numeric value \eqn{>0} may be 
#' specified, with higher values yielding a sparser ensemble. To evaluate the 
#' trade-off between accuracy and sparsity of the final ensemble, inspect
#' \code{pre_object$glmnet.fit} and \code{plot(pre_object$glmnet.fit)}.
#' @param digits Number of decimal places to print
#' @param ... Additional arguments, currently not used.
#' @return Prints information about the fitted prediction rule ensemble.
#' @details Note that the cv error is estimated with data that was also used 
#' for learning rules and may be too optimistic. Use \code{\link{cvpre}} to 
#' obtain a more realistic estimate of future prediction error.
#' @examples \donttest{set.seed(42)
#' airq.ens <- pre(Ozone ~ ., data = airquality[complete.cases(airquality),])
#' print(airq.ens)}
#' @method print pre
#' @seealso \code{\link{pre}}, \code{\link{summary.pre}}, \code{\link{plot.pre}}, 
#' \code{\link{coef.pre}}, \code{\link{importance}}, \code{\link{predict.pre}}, 
#' \code{\link{interact}}, \code{\link{cvpre}} 
#' @export
print.pre <- function(x, penalty.par.val = "lambda.1se", 
                      digits = getOption("digits"),
                      ...) {
  
  if (!(class(x) == "pre" || class(x) == "gpe")) {
    stop("Argument 'x' should be of class 'pre'.")
  }
  
  if (!(length(penalty.par.val) == 1)) {
    stop("Argument 'penalty.par.val' should be a vector of length 1.")
  } else if (!(penalty.par.val == "lambda.min" || 
               penalty.par.val == "lambda.1se" || 
               (is.numeric(penalty.par.val) && penalty.par.val >= 0))) {
    stop("Argument 'penalty.par.val' should be equal to 'lambda.min', 'lambda.1se' or a numeric value >= 0.")
  }
  
  if (!(length(digits) == 1 && digits == as.integer(digits))) {
    stop("Argument 'digits' should be a single integer.")
  }
  
  # function to round values:
  rf <- function(x) {
    signif(x, digits)
  }
  
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
  cat("\n  number of terms = ", x$glmnet.fit$glmnet.fit$df[lambda_ind], 
      "\n  mean cv error (se) = ", rf(x$glmnet.fit$cvm[lambda_ind]), 
      " (", rf(x$glmnet.fit$cvsd[lambda_ind]), ")", "\n\n  cv error type : ",
      x$glmnet.fit$name, "\n\n", sep = "")
  coefs <- coef(x, penalty.par.val = penalty.par.val)
  if (x$family %in% c("gaussian", "poisson", "binomial", "cox")) {
    coefs <- coefs[coefs$coefficient != 0, ]
  } else if (x$family %in% c("mgaussian", "multinomial")) {
    coef_inds <- names(coefs)[!names(coefs) %in% c("rule", "description")]
    coefs <- coefs[rowSums(coefs[,coef_inds]) != 0, ]    
  }
  # always put intercept first:
  is_intercept <- 
    if (is.null(coefs$rule)) {
      rownames(coefs) == "(Intercept)"
    } else {
      coefs$rule == "(Intercept)"
    }
  coefs <- rbind(coefs[is_intercept,], coefs[!is_intercept,])
  
  ## TODO: digits argument appears not to work for multivariate outcomes
  print(coefs, print.gap = 2, quote = FALSE, row.names = FALSE, digits = digits)
  invisible(coefs)
}



#' Summary method for objects of class pre
#'
#' \code{summary.pre} prints information about the generated prediction rule 
#' ensemble to the command line
#' 
#' @param object An object of class \code{\link{pre}}.
#' @inheritParams print.pre
#' @param ... Additional arguments, currently not used.
#' @return Prints information about the fitted prediction rule ensemble.
#' @details Note that the cv error is estimated with data that was also used 
#' for learning rules and may be too optimistic. Use \code{\link{cvpre}} to 
#' obtain a more realistic estimate of future prediction error.
#' @examples \donttest{set.seed(42)
#' airq.ens <- pre(Ozone ~ ., data = airquality[complete.cases(airquality),])
#' summary(airq.ens)}
#' @method summary pre
#' @seealso \code{\link{pre}}, \code{\link{print.pre}}, \code{\link{plot.pre}}, 
#' \code{\link{coef.pre}}, \code{\link{importance}}, \code{\link{predict.pre}}, 
#' \code{\link{interact}}, \code{\link{cvpre}} 
#' @export
summary.pre <- function(object, penalty.par.val = "lambda.1se", ...) {
  
  if (class(object) != "pre") {
    stop("Argument 'object' should be of class 'pre'.")
  }
  
  if (!(length(penalty.par.val) == 1L)) {
    stop("Argument 'penalty.par.val' should be a vector of length 1.")
  } else if (!penalty.par.val %in% c("lambda.min", "lambda.1se") && 
               !(is.numeric(penalty.par.val) && penalty.par.val >= 0)) {
    stop("Argument 'penalty.par.val' should be equal to 'lambda.min', 'lambda.1se' or a numeric value >= 0.")
  }
  
  if (penalty.par.val == "lambda.1se") {
    lambda_ind <- which(object$glmnet.fit$lambda == object$glmnet.fit$lambda.1se)
    cat("\nFinal ensemble with cv error within 1se of minimum: \n  lambda = ", 
        object$glmnet.fit$lambda[lambda_ind])
  }
  if (penalty.par.val == "lambda.min") {
    lambda_ind <- which(object$glmnet.fit$lambda == object$glmnet.fit$lambda.min)
    cat("Final ensemble with minimum cv error: \n\n  lambda = ", 
        object$glmnet.fit$lambda[lambda_ind])
  }
  if (is.numeric(penalty.par.val)) {
    lambda_ind <- which(abs(object$glmnet.fit$lambda - penalty.par.val) == min(abs(
      object$glmnet.fit$lambda - penalty.par.val)))
    cat("Final ensemble with lambda = ", object$glmnet.fit$lambda[lambda_ind])
  }
  cat("\n  number of terms = ", object$glmnet.fit$nzero[lambda_ind], 
      "\n  mean cv error (se) = ", object$glmnet.fit$cvm[lambda_ind], 
      " (", object$glmnet.fit$cvsd[lambda_ind], ")", "\n\n  cv error type : ",
      object$glmnet.fit$name, "\n\n", sep = "")
}


#' Full k-fold cross validation of a prediction rule ensemble (pre)
#' 
#' \code{cvpre} performs k-fold cross validation on the dataset used to create 
#' the specified prediction rule ensemble, providing an estimate of predictive 
#' accuracy on future observations.
#' 
#' @param object An object of class \code{\link{pre}}.
#' @param k integer. The number of cross validation folds to be used.
#' @param verbose logical. Should progress of the cross validation be printed 
#' to the command line?
#' @inheritParams print.pre
#' @param pclass numeric. Only used for binary classification. Cut-off value for the 
#' predicted probabilities that should be used to classify observations to the
#' second class. 
#' @param foldids numeric vector of \code{length(nrow(object$data))} (the number of
#' observations in the training data used to fit the original ensemble). Defaults to
#' \code{NULL}, resulting in the original training observations being randomly 
#' assigned to one of the \eqn{k} folds. Depending on sample size, the number of 
#' factors in the data, the number of factor levels and their distributions, the
#' default may yield errors. See 'Details'. 
#' @param parallel logical. Should parallel foreach be used? Must register parallel 
#' beforehand, such as doMC or others.
#' @param print logical. Should accuracy estimates be printed to the command line?
#' @return Calculates cross-validated estimates of predictive accuracy and prints 
#' these to the command line. For survival regression, accuracy is not calculated, 
#' as there is currently no agreed-upon way to best quantify accuracy in survival 
#' regression models. Users can compute their own accuracy estimates using the 
#' (invisibly returned) cross-validated predictions (\code{$cvpreds}). 
#' Invisibly, a list of three objects is returned: 
#' \code{accuracy} (containing accuracy estimates), \code{cvpreds}
#' (containing cross-validated predictions) and \code{fold_indicators} (a vector indicating
#' the cross validation fold each observation was part of). For (multivariate) continuous 
#' outcomes, accuracy is a list with elements \code{$MSE} (mean squared error on test 
#' observations) and \code{$MAE} (mean absolute error on test observations). For 
#' (binary and multiclass) classification, accuracy is a list with elements 
#' \code{$SEL} (mean squared error on predicted probabilities), \code{$AEL} (mean absolute 
#' error on predicted probabilities), \code{$MCR} (average misclassification error rate) 
#' and \code{$table} (proportion table with (mis)classification rates).
#' @details The random sampling employed by default may yield folds including all 
#' observations with a given level of a given factor. This results in an error, 
#' as it requires predictions for factor levels to be computed that were not 
#' observed in the training data, which is impossible. By manually specifying the
#' \code{foldids} argument, users can make sure all class levels are represented in
#' each of the \eqn{k} training partitions. 
#' @examples \donttest{set.seed(42)
#' airq.ens <- pre(Ozone ~ ., data = airquality[complete.cases(airquality),])
#' airq.cv <- cvpre(airq.ens)}
#' @seealso \code{\link{pre}}, \code{\link{plot.pre}}, 
#' \code{\link{coef.pre}}, \code{\link{importance}}, \code{\link{predict.pre}}, 
#' \code{\link{interact}}, \code{\link{print.pre}} 
#' @export
cvpre <- function(object, k = 10, penalty.par.val = "lambda.1se", pclass = .5, 
                  foldids = NULL, verbose = FALSE, parallel = FALSE,
                  print = TRUE) {
  
  ## check if proper object argument is specified:
  if (class(object) != "pre") {
    stop("Argument 'object' should supply an object of class 'pre'")
  }
  
  ## Check if proper k argument is specified:  
  if (!(length(k) == 1L && k == as.integer(k))) {
    stop("Argument 'k' should be a single positive integer.")
  }  
  
  ## Check if proper verbose argument is specified:  
  if (!(is.logical(verbose) && length(verbose) == 1L)) {
    stop("Argument 'verbose' should be TRUE or FALSE.")
  }  
  
  ## check if pclass is a numeric vector of length 1, and <= 1 and > 0
  if (!(is.numeric(pclass) && length(pclass) == 1L && pclass <= 1 && pclass > 0)) {
    stop("Argument 'verbose' should be TRUE or FALSE.")
  }  
  
  ## check if proper penalty.par.val argument is specified:
  if (!(length(penalty.par.val) == 1L)) {
    stop("Argument 'penalty.par.val' should be a vector of length 1.")
  } else if (!(penalty.par.val == "lambda.min" || 
               penalty.par.val == "lambda.1se" || 
               (is.numeric(penalty.par.val) && penalty.par.val >= 0))) {
    stop("Argument 'penalty.par.val' should be equal to 'lambda.min', 'lambda.1se' or a numeric value >= 0")
  }
  
  ## check if proper parallel argument is specified:
  if (!(is.logical(parallel) && length(parallel) == 1L)) {
    stop("Argument 'parallel' should be TRUE or FALSE")
  }

  ## check is proper foldids argument is specified:
  if (!is.null(foldids)) {
    if (length(foldids) != nrow(object$data)) {
      stop("Argument 'foldids' has length ", length(foldids), ", but should have length ", nrow(object$data), ".")
    } else if (all.equal(foldids, as.integer(foldids))) {
      stop("Argument 'foldids' should be an integer vector, but is not.")
    } 
  }

  ## Set up fold-ids, seeds and object for collecting CV predictions:
  if (is.null(foldids)) {
    foldids <- sample(rep(1:k, length.out = nrow(object$data)), 
                      size = nrow(object$data), replace = FALSE)
  }
  seeds <- sample(k*99, size = k)
  y_ncol <- ifelse(object$family == "multinomial", 
                   nlevels(object$data[,object$y_names]), 
                   length(object$y_names)) 
  cvpreds <- replicate(n = y_ncol, rep(NA, times = nrow(object$data)))
  cl <- object$call
  cl$verbose <- FALSE
  cl$formula <- object$formula 
  
  ## Perform the CV:
  if (parallel) {
    cvpreds_unsorted <- foreach::foreach(i = 1:k, .packages = "pre") %dopar% {
      cl$data <- object$data[foldids != i,]
      set.seed(seeds[i])
      cvobject <- eval(cl)
      predict(cvobject, type = "response", newdata = object$data[foldids == i,], 
              penalty.par.val = penalty.par.val)
    }
    for (i in 1:k) {
      cvpreds[foldids == i,] <- cvpreds_unsorted[[i]]
    }
  } else {
    if (verbose) {
      cat("Running cross validation in fold ")
    }
    for (i in 1:k) {
      
      if (verbose) {
        cat(i, " of ", k, ", ", sep = "")
      }
      
      cl$data <- object$data[foldids != i,]
      set.seed(seeds[i])
      cvobject <- eval(cl)
      cvpreds[foldids == i, ] <- predict(
        cvobject, newdata = object$data[foldids == i,], 
        type = "response", penalty.par.val = penalty.par.val)
      
      if (verbose && i == k) {
        cat("done!\n")
      }
      
    }
  }
  
  ## Collect results:
  accuracy <- list()
  sqrt_N <- sqrt(length(cvpreds) - sum(is.na(cvpreds)))
  if (object$family == "binomial") {
    observed <- object$data[ , object$y_names]
    y_obs <- as.numeric(observed) - 1
    accuracy$SEL<- c(SEL = mean((y_obs - cvpreds)^2, na.rm = TRUE),
      se = sd((y_obs - cvpreds)^2, na.rm = TRUE) / sqrt_N)
    accuracy$AEL <- c(AEL = mean(abs(y_obs - cvpreds), na.rm = TRUE),
      se = sd(abs(y_obs - cvpreds), na.rm = TRUE) / sqrt_N)
    predicted <- factor(cvpreds > pclass)
    levels(predicted) <- levels(observed) 
    accuracy$MCR <- 1 - sum(diag(prop.table(table(predicted, observed))))
    accuracy$table <- prop.table(table(predicted, observed))
  } else if (object$family %in% c("gaussian", "poisson")) {
    y_obs <- object$data[ , object$y_names]
    accuracy$MSE <- c(MSE = mean((y_obs - cvpreds)^2, na.rm = TRUE),
                      se = sd((y_obs - cvpreds)^2, na.rm = TRUE) / sqrt_N)
    accuracy$MAE <- c(MAE = mean(abs(y_obs - cvpreds), na.rm = TRUE),
                      se = sd(abs(y_obs - cvpreds), na.rm = TRUE) / sqrt_N)
  } else if (object$family == "cox") {
    accuracy <- NULL
  } else if (object$family == "mgaussian") {
    y_obs <- object$data[ , object$y_names]
    colnames(cvpreds) <- object$y_names
    accuracy$MSE <- data.frame(MSE = colMeans((y_obs - cvpreds)^2, na.rm = TRUE),
                               se = apply((y_obs - cvpreds)^2, 2, sd, na.rm = TRUE) / sqrt_N)
    accuracy$MAE <- data.frame(MAE = colMeans(abs(y_obs - cvpreds), na.rm = TRUE),
                               se = apply(abs(y_obs - cvpreds), 2, sd, na.rm = TRUE) / sqrt_N)
  } else if (object$family == "multinomial") {
    observed <- object$data[ , object$y_names]
    colnames(cvpreds) <- levels(observed)
    y_obs <- model.matrix( ~ observed + 0)
    colnames(y_obs) <- levels(observed)
    accuracy$SEL<- data.frame(SEL = colMeans((y_obs - cvpreds)^2, na.rm = TRUE),
                              se = apply((y_obs - cvpreds)^2, 2, sd, na.rm = TRUE) / sqrt_N)
    accuracy$AEL <- data.frame(AEL = colMeans(abs(y_obs - cvpreds), na.rm = TRUE),
                               se = apply(abs(y_obs - cvpreds), 2, sd, na.rm = TRUE) / sqrt_N)
    predicted <- factor(apply(cvpreds, 1, function(x) which(x == max(x))), levels = 1:ncol(cvpreds))
    levels(predicted) <- levels(observed)
    accuracy$MCR <- 1 - sum(diag(prop.table(table(predicted, observed))))
    accuracy$table <- prop.table(table(predicted, observed))
  }
  if (print && object$family != "cox") print(accuracy)
  return(invisible(list(accuracy = accuracy, cvpreds = cvpreds, fold_indicators = foldids)))
}






#' Coefficients for the final prediction rule ensemble
#'
#' \code{coef.pre} returns coefficients for prediction rules and linear terms in 
#' the final ensemble
#' 
#' @param object object of class \code{\link{pre}}
#' @inheritParams print.pre
#' @param ... additional arguments to be passed to \code{\link[glmnet]{coef.glmnet}}.
#' @return returns a dataframe with 3 columns: coefficient, rule (rule or 
#' variable name) and description (\code{NA} for linear terms, conditions for 
#' rules).
#' @details In some cases, duplicated variable names may appear in the model.
#' For example, the first variable is a factor named 'V1' and there are also
#' variables named 'V10' and/or 'V11' and/or 'V12' (etc). Then for 
#` selecting the final ensemble, if linear terms are also included,
#' for the binary factor V1, dummy contrast variables will be created, named 
#' 'V10', 'V11', 'V12' (etc). As should be clear from this example, this yields 
#' duplicated variable names, which may yield problems, for example in the 
#' calculation of predictions and importances, later on. This can be prevented 
#' by renaming factor variables with numbers in their name, prior to analysis.
#' 
#' 
#' @examples \donttest{set.seed(42)
#' airq.ens <- pre(Ozone ~ ., data = airquality[complete.cases(airquality),])
#' coefs <- coef(airq.ens)}
#' @method coef pre
#' @seealso \code{\link{pre}}, \code{\link{plot.pre}}, 
#' \code{\link{cvpre}}, \code{\link{importance}}, \code{\link{predict.pre}}, 
#' \code{\link{interact}}, \code{\link{print.pre}} 
#' @export
coef.pre <- function(object, penalty.par.val = "lambda.1se", ...)
{
  
  ## TODO: Add argument on whether learners with zero coefficients should
  ## be included ot not
  
  ## check if proper object argument is specified:
  if(class(object) != "pre") {
    stop("Argument 'object' should supply an object of class 'pre'")
  }

  ## check if proper penalty.par.val argument is specified:
  if (!(length(penalty.par.val) == 1)) {
    stop("Argument 'penalty.par.val' should be a vector of length 1.")
  } else if (!(penalty.par.val == "lambda.min" || 
               penalty.par.val == "lambda.1se" || 
               (is.numeric(penalty.par.val) && penalty.par.val >= 0))) {
    stop("Argument 'penalty.par.val' should be equal to 'lambda.min', 'lambda.1se' or a numeric value >= 0")
  }
  
  if (object$family %in% c("gaussian", "binomial", "poisson", "cox")) {
    coefs <- as(coef.glmnet(object$glmnet.fit, s = penalty.par.val, ...), 
                Class = "matrix")
  } else if (object$family %in% c("mgaussian", "multinomial")) {
    coefs <- sapply(coef(object$glmnet.fit, s = penalty.par.val, ...), as, 
                    Class = "matrix")
    rownames(coefs) <- rownames(coef(object$glmnet.fit)[[1]])
  }
  
  rownames(coefs) <- gsub("`", "", rownames(coefs))
  
  # coefficients for normalized variables should be unnormalized: 
  if (object$normalize & !is.null(object$x_scales) & object$type != "rules") {
    coefs[names(object$x_scales),] <- coefs[names(object$x_scales),] /
      object$x_scales
  }
  if (object$family %in% c("gaussian", "binomial", "poisson", "cox")) {
    coefs <- data.frame(coefficient = coefs[,1], rule = rownames(coefs), 
                        stringsAsFactors = FALSE)
  } else if (object$family %in% c("mgaussian", "multinomial")) {
    coefs <- data.frame(coefficient = coefs, rule = rownames(coefs), 
                        stringsAsFactors = FALSE)
  }
  
  # check whether there's duplicates in the variable names:
  # (can happen, for example, due to labeling of dummy indicators for factors)
  if (!(length(unique(coefs$rule)) == length(coefs$rule))) { 
    replicates_in_variable_names <- TRUE
    warning("There are variables in the model with overlapping variable names. This may result in errors, or results may not be valid. If predictor variables of type factor were specified with numbers in their name, consider renaming these. See 'Details' under ?coef.pre.") 
  } else {
    replicates_in_variable_names <- FALSE
  }
  if (object$type != "linear" && !is.null(object$rules)) {
    # We set sort to FALSE to get comparable results across platforms
    coefs <- base::merge.data.frame(coefs, object$rules, all.x = TRUE, sort = FALSE)
    coefs$description <- as.character(coefs$description)
  } else {
    if (object$family %in% c("mgaussian", "multinomial")) {
      coefs <- data.frame(rule = coefs$rule, 
                          description = rep(NA, times = nrow(coefs)), 
                          coefs[,which(names(coefs) != "rule")],
                          stringsAsFactors = FALSE)      
    } else {
      coefs <- data.frame(rule = coefs$rule, 
                          description = rep(NA, times = nrow(coefs)), 
                          coefficient = coefs[,1],
                          stringsAsFactors = FALSE)
    }
  }
  
  ## Description of the intercept should be 1:
  coefs$description[which(coefs$rule == "(Intercept)")] <- "1"
  
  ## Description of input variables:
  coefs$description[is.na(coefs$description)] <- coefs$rule[is.na(coefs$description)]

  # include winsorizing points in the description if they were used in 
  # generating the ensemble (and if there are no duplicate variable names):  
  if (!is.null(object$wins_points) && !replicates_in_variable_names) { 
    wp <- object$wins_points[!is.na(object$wins_points$value), ]
    coefs[coefs$rule %in% wp$varname, ][
      order(coefs[coefs$rule %in% wp$varname,]$rule), ]$description <- 
      wp[order(wp$varname), ]$value  
  }
  
  if (object$family %in% c("gaussian", "binomial", "poisson", "cox")) {
    return(coefs[order(abs(coefs$coefficient), decreasing = TRUE),])
  } else if (object$family %in% c("mgaussian", "multinomial")) {
    return(coefs[order(abs(coefs[,3]), decreasing = TRUE),])    
  }
}




#' Predicted values based on final prediction rule ensemble
#'
#' \code{predict.pre} generates predictions based on the final prediction rule
#' ensemble, for training or new (test) observations
#'
#' @param object object of class \code{\link{pre}}.
#' @param newdata optional dataframe of new (test) observations, including all
#' predictor variables used for deriving the prediction rule ensemble.
#' @inheritParams print.pre
#' @param type character string. The type of prediction required; the default
#' \code{type = "link"} is on the scale of the linear predictors. Alternatively,
#' for count and factor outputs, \code{type = "response"} may be specified to obtain
#' the fitted mean and fitted probabilities, respectively; \code{type = "class"} 
#' returns the predicted class membership.
#' @param ... further arguments to be passed to 
#' \code{\link[glmnet]{predict.cv.glmnet}}.
#' @details If \code{newdata} is not provided, predictions for training data will be 
#' returned.
#' @examples \donttest{set.seed(1)
#' train <- sample(1:sum(complete.cases(airquality)), size = 100)
#' set.seed(42)
#' airq.ens <- pre(Ozone ~ ., data = airquality[complete.cases(airquality),][train,])
#' predict(airq.ens)
#' predict(airq.ens, newdata = airquality[complete.cases(airquality),][-train,])}
#' @import Matrix
#' @method predict pre
#' @seealso \code{\link{pre}}, \code{\link{plot.pre}}, 
#' \code{\link{coef.pre}}, \code{\link{importance}}, \code{\link{cvpre}}, 
#' \code{\link{interact}}, \code{\link{print.pre}}, 
#' \code{\link[glmnet]{predict.cv.glmnet}}
#' @export
predict.pre <- function(object, newdata = NULL, type = "link",
                        penalty.par.val = "lambda.1se", ...)
{
  ## Check if proper object argument is specified:
  if(class(object) != "pre") {
    stop("Argument 'object' should supply an object of class 'pre'")
  }
  
  ## check if proper type argument is specified:
  if (length(type) != 1L || !is.character(type)) {
    stop("Argument 'type' should be a character vector of length 1")
  }
  
  ## check if proper penalty.par.val argument is specified:
  if (!(length(penalty.par.val) == 1L)) {
    stop("Argument 'penalty.par.val' should be a vector of length 1.")
  } else if (!(penalty.par.val%in% c("lambda.min", "lambda.1se")) && 
               !(is.numeric(penalty.par.val) && penalty.par.val >= 0)) {
    stop("Argument 'penalty.par.val' should be equal to 'lambda.min', 'lambda.1se' or a numeric value >= 0")
  }

  if (is.null(newdata)) {
    
    newdata <- object$modmat
  
  } else {

    ## Have to prepare newdata for get_modmat():
    
    ## check if proper newdata argument is specified, if specified:    
    if (!is.data.frame(newdata)) {
      stop("newdata should be a data frame.")
    }
    
    ## Get winsfrac (to pass on to get_modmat later):
    winsfrac <- (object$call)$winsfrac
    if(is.null(winsfrac))
      winsfrac <- formals(pre)$winsfrac
    
    ## Check if variable names and classes are the same in newdata as in object$data:
    if (!all(object$x_names %in% names(newdata))) {
      newdata <- model.frame(as.Formula((object$call)$formula), data = newdata, 
                             rhs = NULL, lhs = 0, na.action = NULL)
    } else {
      newdata <- newdata[ , object$x_names]
    }
    
    ## Coerce character and logical variables to factors:
    if (any(char_names <- sapply(newdata, is.character))) {
      char_names <- names(newdata)[char_names]
      data[ , char_names] <- sapply(newdata[ , char_names], factor)
    }
    if (any(logic_names <- sapply(newdata, is.logical))) {
      logic_names <- names(newdata)[logic_names]
      newdata[ , logic_names] <- sapply(newdata[ , logic_names], factor)
    } 
    
    ## Coerce ordered categorical variables to numeric, if necessary:
    if (if (is.null((object$call)$ordinal)) {
      formals(pre)$ordinal
    } else {
      (object$call)$ordinal
    }) {
      if (any(ord_var_inds <- sapply(newdata, is.ordered))) {
        newdata[ , ord_var_inds] <- sapply(newdata[ , ord_var_inds], as.numeric)
      }
    }
    
    if (any(is.na(newdata))) {
      newdata <- newdata[complete.cases(newdata),]
      warning("Some observations in newdata have missing predictor variable values and will be removed.", immediate. = TRUE)
    }
      
    ## Check and set factor levels of newdata to variable levels in object$data:
    if (any(factor_inds <- sapply(newdata, is.factor))) {
      for (i in names(newdata)[factor_inds]) {
        if (all(levels(newdata[ , i]) %in% levels(object$data[ , i]))) {
          levels(newdata[ , i]) <- levels(object$data[ , i])
        } else {
          stop("Variable ", i, " has levels not present in training data. Cannot compute predictions.")
        }
      }
    }
    
    newdata <- get_modmat(
      wins_points = object$wins_points, 
      x_scales = object$x_scales, 
      formula = object$formula, 
      data = newdata, 
      rules = if (object$type == "linear" || is.null(object$rules)) {NULL} else {
        structure(object$rules$description, names = object$rules$rule)}, 
      type = object$type, 
      winsfrac = winsfrac,
      x_names = object$x_names, 
      normalize = object$normalize,
      y_names = NULL,
      confirmatory = object$call$confirmatory)
    
    newdata <- newdata$x
  }
  
  ## Get predictions:
  if (object$family %in% c("gaussian", "binomial", "poisson", "cox")) {
    preds <- predict.cv.glmnet(object$glmnet.fit, newx = newdata, 
                               s = penalty.par.val, type = type, ...)[,1]
  } else if (object$family %in% c("mgaussian", "multinomial")) {
    if (object$family == "multinomial" && type == "class") {
      preds <- predict.cv.glmnet(object$glmnet.fit, newx = newdata, 
                                 s = penalty.par.val, type = type, ...)[,1]
    } else {
      preds <- predict.cv.glmnet(object$glmnet.fit, newx = newdata, 
                                 s = penalty.par.val, type = type, ...)[,,1]
    }

  }
  return(preds)
}




#' Create partial dependence plot for a single variable in a prediction rule 
#' ensemble (pre)
#'
#' \code{singleplot} creates a partial dependence plot, which shows the effect of
#' a predictor variable on the ensemble's predictions. Note that plotting partial 
#' dependence is computationally intensive. Computation time will increase fast 
#' with increasing numbers of observations and variables. For large 
#' datasets, package `plotmo` (Milborrow, 2019) provides more efficient functions 
#' for plotting partial dependence and also supports `pre` models. 
#'
#' @param object an object of class \code{\link{pre}}
#' @param varname character vector of length one, specifying the variable for
#' which the partial dependence plot should be created. Note that \code{varname}
#' should correspond to the variable as described in the model formula used
#' to generate the ensemble (i.e., including functions applied to the variable).
#' @inheritParams print.pre
#' @param nvals optional numeric vector of length one. For how many values of x
#' should the partial dependence plot be created?
#' @param type character string. Type of prediction to be plotted on y-axis.
#' \code{type = "response"} gives fitted values for continuous outputs and
#' fitted probabilities for nominal outputs. \code{type = "link"} gives fitted
#' values for continuous outputs and linear predictor values for nominal outputs.
#' @param ylab character. Label to be printed on the y-axis.
#' @param ... Further arguments to be passed to 
#' \code{\link[graphics]{plot.default}}.
#' @details By default, a partial dependence plot will be created for each unique
#' observed value of the specified predictor variable. When the number of unique
#' observed values is large, this may take a long time to compute. In that case,
#' specifying the \code{nvals} argument can substantially reduce computing time. When the
#' \code{nvals} argument is supplied, values for the minimum, maximum, and \code{(nvals - 2)}
#' intermediate values of the predictor variable will be plotted. Note that \code{nvals}
#' can be specified only for numeric and ordered input variables. If the plot is
#' requested for a nominal input variable, the \code{nvals} argument will be
#' ignored and a warning printed.
#' 
#' See also section 8.1 of Friedman & Popescu (2008).
#' 
#' @references Fokkema, M. (2018). Fitting prediction rule ensembles with R package pre.
#' \url{https://arxiv.org/abs/1707.07149}.
#' 
#' Friedman, J. H., & Popescu, B. E. (2008). Predictive learning 
#' via rule ensembles. \emph{The Annals of Applied Statistics, 2}(3), 916-954.
#' 
#' Milborrow, S. (2019). plotmo: Plot a model's residuals, response, and partial 
#' dependence plots. \url{https://CRAN.R-project.org/package=plotmo}
#' 
#' @examples \donttest{set.seed(42)
#' airq.ens <- pre(Ozone ~ ., data = airquality[complete.cases(airquality),])
#' singleplot(airq.ens, "Temp")}
#' @seealso \code{\link{pre}}, \code{\link{pairplot}}
#' @export
singleplot <- function(object, varname, penalty.par.val = "lambda.1se",
                       nvals = NULL, type = "response", ylab = "predicted", 
                       ...)
{
 
  ## Check family:
  if (object$family %in% c("mgaussian", "multinomial")) {
    stop("Function singleplot not implemented yet for multivariate and multinomial outcomes.")
  }
  
  ## Check if proper object argument is specified: 
  if(class(object) != "pre") {
    stop("Argument 'object' should be an object of class 'pre'")
  }
  
  ## Check if proper varname argument is specified: 
  if (length(varname) != 1L || !is.character(varname)) {
    stop("Argument 'varname' should be a character vector of length 1.")
  } else if (!(varname %in% object$x_names)) {
    varnames <- grep(varname, x = object$x_names, value = TRUE, fixed = TRUE)
    if (length(varnames > 0)) {
      stop("Argument 'varname' should specify the variable name as specified in the model formula (e.g., ", paste0(paste0("'", varnames, "'"), collapse = " or "), ").")
    } else {
      stop("Argument 'varname' should specify the name of a variable used to generate the ensemble.")
    }
  }
  
  ## Check if proper penalty.par.val argument is specified: 
  if (length(penalty.par.val) != 1L) {
    stop("Argument 'penalty.par.val' should be a vector of length 1.")
  } else if (!( penalty.par.val %in% c("lambda.min", "lambda.1se") || 
               (is.numeric(penalty.par.val) && penalty.par.val >= 0))) {
    stop("Argument 'penalty.par.val' should be equal to 'lambda.min', 'lambda.1se' or a numeric value >= 0")
  }
  
  ## Check if proper nvals argument is specified: 
  if (!is.null(nvals)) {
    if (length(nvals) != 1L || nvals != as.integer(nvals)) {
      stop("Argument 'nvals' should be an integer vector of length 1.")
    } else if (is.factor(object$data[,varname]) && !is.null(nvals)) {
      warning("Plot is requested for variable of class factor. Value specified for
              nvals will be ignored.", immediate. = TRUE)
      nvals <- NULL
    }
  }
  
  ## Check if proper type argument is specified: 
  if (length(type) != 1L || !is.character(type)) {
    stop("Argument 'type' should be a single character string.")
  }
  
  # Generate expanded dataset:
  if (is.null(nvals)) {
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
    type = "l", ylab = ylab, xlab = varname, ...)
}





#' Create partial dependence plot for a pair of predictor variables in a prediction 
#' rule ensemble (pre)
#'
#' \code{pairplot} creates a partial dependence plot to assess the effects of a
#' pair of predictor variables on the predictions of the ensemble. Note that plotting 
#' partial dependence is computationally intensive. Computation time will increase 
#' fast with increasing numbers of observations and variables. For large 
#' datasets, package `plotmo` (Milborrow, 2019) provides more efficient functions 
#' for plotting partial dependence and also supports `pre` models. 
#'
#' @param object an object of class \code{\link{pre}}
#' @param varnames character vector of length two. Currently, pairplots can only
#' be requested for non-nominal variables. If varnames specifies the name(s) of
#' variables of class \code{"factor"}, an error will be printed.
#' @inheritParams print.pre
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
#' See also section 8.1 of Friedman & Popescu (2008).
#' 
#' @note Function \code{pairplot} uses package akima to construct interpolated 
#' surfaces and  has an ACM license that restricts applications to non-commercial 
#' usage, see 
#' \url{https://www.acm.org/publications/policies/software-copyright-notice}
#' Function \code{pairplot} prints a note referring to this ACM licence.
#' @examples \donttest{set.seed(42)
#' airq.ens <- pre(Ozone ~ ., data = airquality[complete.cases(airquality),])
#' pairplot(airq.ens, c("Temp", "Wind"))}
#' @export
#' @references Fokkema, M. (2018). Fitting prediction rule ensembles with R package pre.
#' \url{https://arxiv.org/abs/1707.07149}.
#' 
#' Friedman, J. H., & Popescu, B. E. (2008). Predictive learning 
#' via rule ensembles. \emph{The Annals of Applied Statistics, 2}(3), 916-954.
#' 
#' Milborrow, S. (2019). plotmo: Plot a model's residuals, response, and partial 
#' dependence plots. \url{https://CRAN.R-project.org/package=plotmo}
#' @import graphics
#' @seealso \code{\link{pre}}, \code{\link{singleplot}} 
#' @export
pairplot <- function(object, varnames, type = "both", 
                     penalty.par.val = "lambda.1se", 
                     nvals = c(20, 20), pred.type = "response", ...)
{
  
  ## Check if proper object argument is specified: 
  if(class(object) != "pre") {
    stop("Argument 'object' should be an object of class 'pre'")
  }
  
  ## Check family:
  if (object$family %in% c("mgaussian", "multinomial")) {
    stop("Function pairplot not implemented yet for multivariate and multinomial outcomes.")
  }
  
  ## Check if proper varnames argument is specified: 
  if (length(varnames) != 2L || !is.character(varnames)) {
    stop("Argument 'varnames' should be a character vector of length 2.")
  } else if (!(all(varnames %in% object$x_names))) {
        varname <- grep(varnames[1], x = object$x_names, value = TRUE, fixed = TRUE)
        varnames <- c(varname, grep(varnames[2], x = object$x_names, value = TRUE, fixed = TRUE))
    if (length(varnames > 0)) {
      stop("Argument 'varnames' should specify the variable names as specified in the model formula (e.g., ", paste0(paste0("'", varnames, "'"), collapse = " and/or "), ").")
    } else {
      stop("Argument 'varnames' should specify names of variables used to generate the ensemble.")
    }
  } else if (any(sapply(object$data[,varnames], is.factor))) {
    stop("3D partial dependence plots are currently not supported for factors.")
  }
  
  ## Check if proper penalty.par.val argument is specified: 
  if (!(length(penalty.par.val) == 1)) {
    stop("Argument 'penalty.par.val' should be a vector of length 1.")
  } else if (!(penalty.par.val %in% c("lambda.min", "lambda.1se") || 
               (is.numeric(penalty.par.val) && penalty.par.val >= 0))) {
    stop("Argument 'penalty.par.val' should be equal to 'lambda.min', 'lambda.1se' or a numeric value >= 0")
  }
  
  ## Check if proper nvals argument is specified: 
  if (!(length(nvals) == 2 && nvals == as.integer(nvals))) {
    stop("Argument 'nvals' should be an integer vector of length 2.")
  }
  
  ## Check if proper type argument is specified: 
  if (!(length(type) == 1 && is.character(type))) {
    stop("Argument 'type' should be equal to 'heatmap', 'contour', 'both' or 'perspective'.")
  }
  
  ## Check if proper pred.type argument is specified: 
  if (!(length(type) == 1 && is.character(type))) {
    stop("Argument 'type' should be a single character string.")
  }
  
  ## check if pakcage akima is installed:
  if (!(requireNamespace("akima"))) {
    stop("Function pairplot requires package akima. Download and install package
         akima from CRAN, and run again.")
  }


  # generate expanded dataset:
  if (is.null(nvals)){
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
  pred_vals <- predict.pre(object, newdata = exp_dataset, type = pred.type,
                           penalty.par.val = penalty.par.val)
  
  # create plot:
  if (is.null(nvals)) {nvals <- 3}
  xyz <- akima::interp(exp_dataset[,varnames[1]], exp_dataset[,varnames[2]],
                       pred_vals, duplicate = "mean")
  if (type == "heatmap" || type == "both") {
    if (is.null(match.call()$col)) {
      colors <- rev(c("#D33F6A", "#D95260", "#DE6355", "#E27449", "#E6833D", 
               "#E89331", "#E9A229", "#EAB12A", "#E9C037", "#E7CE4C", 
               "#E4DC68", "#E2E6BD"))
      image(xyz, xlab = varnames[1], ylab = varnames[2], 
            col = colors, ...)
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


#' Calculate importances of baselearners and input variables in a prediction 
#' rule ensemble (pre)
#'
#' \code{importance} calculates importances for rules, linear terms and input
#' variables in the prediction rule ensemble (pre), and creates a bar plot 
#' of variable importances.
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
#' @inheritParams print.pre
#' @param round integer. Number of decimal places to round numeric results to.
#' If \code{NA} (default), no rounding is performed.
#' @param plot logical. Should variable importances be plotted?
#' @param ylab character string. Plotting label for y-axis. Only used when
#' \code{plot = TRUE}.
#' @param main character string. Main title of the plot. Only used when
#' \code{plot = TRUE}.
#' @param diag.xlab logical. Should variable names be printed diagonally (that
#' is, in a 45 degree angle)? Alternatively, variable names may be printed 
#' vertically by specifying \code{diag.xlab = FALSE} and \code{las = 2}.
#' @param abbreviate integer or logical. Number of characters to abbreviate 
#' x axis names to. If \code{FALSE}, no abbreviation is performed. 
#' @param diag.xlab.hor numeric. Horizontal adjustment for lining up variable
#' names with bars in the plot if variable names are printed diagonally.
#' @param diag.xlab.vert positive integer. Vertical adjustment for position
#' of variable names, if printed diagonally. Corresponds to the number of 
#' character spaces added after variable names. 
#' @param cex.axis numeric. The magnification to be used for axis annotation
#' relative to the current setting of \code{cex}.
#' @param legend logical or character. Should legend be plotted for multinomial
#' or multivariate responses and if so, where? Defaults to \code{"topright"}, 
#' which puts the legend in the top-right corner of the plot. Alternatively, 
#' \code{"bottomright"}, \code{"bottom"}, \code{"bottomleft"}, \code{"left"}, 
#' \code{"topleft"}, \code{"top"}, \code{"topright"}, \code{"right"}, 
#' \code{"center"} and \code{FALSE} (which omits the legend) can be specified.
#' @param ... further arguments to be passed to \code{barplot} (only used
#' when \code{plot = TRUE}).
#' @return A list with two dataframes: \code{$baseimps}, giving the importances 
#' for baselearners in the ensemble, and \code{$varimps}, giving the importances 
#' for all predictor variables.
#' @details See also sections 6 and 7 of Friedman & Popecus (2008).
#' @examples \donttest{set.seed(42)
#' airq.ens <- pre(Ozone ~ ., data = airquality[complete.cases(airquality),])
#' # calculate global importances:
#' importance(airq.ens)
#' # calculate local importances (default: over 25% highest predicted values):
#' importance(airq.ens, global = FALSE)
#' # calculate local importances (custom: over 25% lowest predicted values):
#' importance(airq.ens, global = FALSE, quantprobs = c(0, .25))}
#' @references Fokkema, M. (2018). Fitting prediction rule ensembles with R package pre.
#' \url{https://arxiv.org/abs/1707.07149}.
#' 
#' Friedman, J. H., & Popescu, B. E. (2008). Predictive learning 
#' via rule ensembles. \emph{The Annals of Applied Statistics, 2}(3), 916-954.
#' @seealso \code{\link{pre}}
#' @export
importance <- function(object, standardize = FALSE, global = TRUE,
                       quantprobs = c(.75, 1), penalty.par.val = "lambda.1se", 
                       round = NA, plot = TRUE, ylab = "Importance",
                       main = "Variable importances", abbreviate = 10L,
                       diag.xlab = TRUE, diag.xlab.hor = 0, diag.xlab.vert = 2,
                       cex.axis = 1, legend = "topright", ...)
{

  if (!inherits(object, what = "pre")) {
    stop("Specified object is not of class 'pre'.")
  }
  
  if (!global) {
    if (object$family %in% c("mgaussian", "multinomial")) {
      warning("Local importances cannot be calculated for multivariate and multinomial outcomes. Global importances will be returned.")
      global <- TRUE 
    }
  }
  
  if (standardize && object$family %in% c("multinomial", "binomial", "cox")) {
    warning("Standardized importances cannot be calculated for binary, multinomial or survival responses. Unstandardized importances will be returned.")
    standardize <- FALSE
  }
  
  ## Step 1: Calculate the importances of the base learners:
  
  # get base learner coefficients:
  coefs <- coef(object, penalty.par.val = penalty.par.val)
  if (object$family %in% c("mgaussian", "multinomial")) {
    coef_inds <- names(coefs)[!names(coefs) %in% c("rule", "description")]
  }

  # only continue when there are nonzero terms besides intercept:
  if ((object$family %in% c("gaussian", "binomial", "poisson") && 
      sum(coefs$coefficient != 0) > 1L ) || 
      (object$family %in% c("mgaussian", "multinomial") && 
       sum(rowSums(coefs[,coef_inds]) != 0) > 1L) ||
      (object$family == "cox" && sum(coefs$coefficient != 0) > 0)) { 
    # give factors a description:
    if (any(is.na(coefs$description))) {
      coefs$description[is.na(coefs$description)] <-
        paste0(as.character(coefs$rule)[is.na(coefs$description)], " ")
    }
    coefs <- coefs[order(coefs$rule),]
    # Get sds for every baselearner:
    if (global) {
      # object$x_scales should be used to get correct SDs for linear terms:
      if (object$family == "cox") {
        sds <- apply(object$modmat, 2, sd, na.rm = TRUE)  
      } else {
        sds <- c(0, apply(object$modmat, 2, sd, na.rm = TRUE))          
      }
      if (standardize) {
        if (object$family == "mgaussian") {
          sd_y <- sapply(object$data[,object$y_names], sd)
        } else if (object$family %in% c("gaussian", "poisson")) { 
          sd_y <- sd(as.numeric(object$data[,object$y_names]))
        }
      }
    } else {
      preds <- predict.pre(object, newdata = object$data, type = "response",
                           penalty.par.val = penalty.par.val)
      local_modmat <- object$modmat[preds >= quantile(preds, probs = quantprobs[1]) &
                                      preds <= quantile(preds, probs = quantprobs[2]),]
      if (nrow(local_modmat) < 2) {stop("Selected subregion contains less than 2 observations, importances cannot be calculated")}
      # object$x_scales should be used to get correct SDs for linear terms:
      if (object$family == "cox") {
        ## cox prop haz model has no intercept, so should be omitted
        sds <- apply(local_modmat, 2, sd, na.rm = TRUE)
      } else {
        sds <- c(0, apply(local_modmat, 2, sd, na.rm = TRUE))
      }

      if (standardize) {
        sd_y <- sd(object$data[preds >= quantile(preds, probs = quantprobs[1]) & 
                                 preds <= quantile(preds, probs = quantprobs[2]),
                               object$y_names])
      }
    }
    
    ## Check if there are any " ` " marks in sd names, if so remove:
    if (any(grepl("`", names(sds), fixed = TRUE))) {
      names(sds) <- gsub("`", "", names(sds), fixed = TRUE)
    }
    
    if(object$normalize) {
      sds[names(object$x_scales)] <- sds[names(object$x_scales)] * object$x_scales
    }
    
    if (object$family != "cox") {
      names(sds)[1] <- "(Intercept)"
    }
    
    sds <- sds[order(names(sds))]
  
    if (any(names(sds) != coefs$rule)) {
      warning("There seems to be a problem with the ordering or size of the coefficient and sd vectors. Importances cannot be calculated.")
    }
    
    # baselearner importance is given by abs(coef*SD) (F&P section 6):
    if (object$family %in% c("multinomial", "mgaussian")) {
      baseimps <- data.frame(coefs, sd = sds)
      baseimps[,gsub("coefficient", "importance", coef_inds)] <- abs(sapply(baseimps[,coef_inds], function(x) {x*sds}))
    } else {
      baseimps <- data.frame(coefs, sd = sds, imp = abs(coefs$coefficient)*sds)
    }
    
    if (standardize) {
      if (object$family == "mgaussian") {
        for (i in gsub("coefficient", "importance", coef_inds)) {
          baseimps[,i] <- baseimps[,i] / sd_y[gsub("importance.", "", i)]
        }
      } else if (object$family %in% c("gaussian", "poisson")) {
        baseimps$imp <- baseimps$imp / sd_y
      }
    }

    ## Remove nonzero terms:
    if (object$family %in% c("mgaussian", "multinomial")) {
      baseimps <- baseimps[rowSums(baseimps[,coef_inds]) != 0, ]   
    } else {
      baseimps <- baseimps[baseimps$coefficient != 0,]
    }
    
    ## Omit intercept:
    baseimps <- baseimps[baseimps$description != "1",]
    
    # Calculate the number of terms in each rule:
    baseimps$nterms <- NA
    for(i in 1:nrow(baseimps)) {
      # If there is " & " in description, there are at least 2 terms/variables 
      # in the base learner:
      if (grepl(" & ", baseimps$description[i])) {
        baseimps$nterms[i] <- length(gregexpr("&", baseimps$description)[[i]]) + 1L
      } else {
        baseimps$nterms[i] <- 1L # if not, the number of terms = 1
      }
    }
    
    
    ## Step 2: Calculate variable importances:
    
    if (object$family %in% c("mgaussian", "multinomial")) {
      varimps <- data.frame(varname = object$x_names, stringsAsFactors = FALSE)
      varimps[,gsub("coefficient", "importance", coef_inds)] <- 0
    } else {
      varimps <- data.frame(varname = object$x_names, imp = 0,
                            stringsAsFactors = FALSE)
    }
    
    for(i in 1:nrow(varimps)) { # for every variable:
      
      ## Get imps from rules and linear terms:
      for(j in 1:nrow(baseimps)) {
        # if the variable name appears in the rule:
        #   (Note: EXACT matches are needed, so 1) there should be a space before 
        #     and after the variable name in the rule and thus 2) there should be 
        #     a space added before the description of the rule)
        if(grepl(paste0(" ", varimps$varname[i], " "), paste0(" ", baseimps$description[j]))) {
          # count the number of times it appears in the rule:
          n_occ <- length(gregexpr(paste0(" ", varimps$varname[i], " "),
                                   paste0(" ", baseimps$description[j]), fixed = TRUE)[[1]])
          # add it to the importance of the variable:
          if (object$family %in% c("mgaussian", "multinomial")) {
            varimps[i, gsub("coefficient", "importance", coef_inds)] <- 
              varimps[i, gsub("coefficient", "importance", coef_inds)] + 
              (n_occ * baseimps[j, gsub("coefficient", "importance", coef_inds)] / baseimps$nterms[j])
          } else {
            varimps$imp[i] <- varimps$imp[i] + (n_occ * baseimps$imp[j] / baseimps$nterms[j])
          }
        }
      }
        
      ## Get imps from factors:
      if (is.factor(object$data[ , varimps$varname[i]])) { # && 
          # !is.ordered(object$data[ , varimps$varname[i]])) {
        ## Sum those baseimps$imp for which baseimps$rule has varimps$varname[i] as part of its name
        if (object$family %in% c("mgaussian", "multinomial")) {
          varimps[i, gsub("coefficient", "importance", coef_inds)] <-
            varimps[i, gsub("coefficient", "importance", coef_inds)] +
            colSums(baseimps[grepl(varimps$varname[i], baseimps$rule, fixed = TRUE), 
                             gsub("coefficient", "importance", coef_inds)])
        } else {
          varimps$imp[i] <- varimps$imp[i] + 
            sum(baseimps$imp[grepl(varimps$varname[i], baseimps$rule, fixed = TRUE)])
        }
      }
    }
    
    
    ## Step 3: Return (and plot) importances:
    

    
    if (object$family %in% c("mgaussian", "multinomial")) {
      varimps <- varimps[rowSums(varimps[ , gsub("coefficient", "importance", coef_inds)]) != 0, ]
      ord <- order(rowSums(varimps[ , gsub("coefficient", "importance", coef_inds)]), 
                   decreasing = TRUE, method = "radix")
      varimps <- varimps[ord, ]
    } else {
      baseimps <- baseimps[order(baseimps$imp, decreasing = TRUE, method = "radix"), ]
      varimps <- varimps[order(varimps$imp, decreasing = TRUE, method = "radix"), ]
      varimps <- varimps[varimps$imp != 0, ]
    }
    
    if (plot & nrow(varimps) > 0) {
      if (is.character(legend)) {
        args.legend <- list(x = legend)
        legend.text <- TRUE
      } else {
        legend.text <- NULL
        args.legend <- NULL
      }
      if (object$family %in% c("mgaussian", "multinomial")) {
        plot_varimps <- t(varimps[ , gsub("coefficient", "importance" , coef_inds)])
        colnames(plot_varimps) <- abbreviate(varimps$varname, minlength = abbreviate)
        rownames(plot_varimps) <- gsub("coefficient.", "" , coef_inds)
        if (diag.xlab) {
          xlab.pos <- barplot(plot_varimps, beside = TRUE, ylab = ylab, 
                              names.arg = rep("", times = ncol(plot_varimps)), 
                              main = main, cex.axis = cex.axis, 
                              legend.text = legend.text, 
                              args.legend = args.legend, ...)
          xlab.pos <- xlab.pos[nrow(xlab.pos),]
          ## add specified number of trailing spaces to variable names:
          plotnames <- varimps$varname
          if (is.numeric(abbreviate) && abbreviate > 0) {
            plotnames <- abbreviate(plotnames, minlength = abbreviate)
          }
          if (diag.xlab.vert > 0) {
            for (i in 1:diag.xlab.vert) {
              plotnames <- paste0(plotnames, " ")
            }
          }
          text(xlab.pos + diag.xlab.hor, par("usr")[3], srt = 45, adj = 1, xpd = TRUE, 
               labels = plotnames, cex = cex.axis)
        } else {
          barplot(plot_varimps, beside = TRUE, main = main, ylab = ylab, 
                  legend.text = legend.text, args.legend = args.legend,
                  ...)
        }
      } else {
        if (diag.xlab) {
          xlab.pos <- barplot(height = varimps$imp, xlab = "", ylab = ylab, 
                              main = main, cex.axis = cex.axis, ...)
          ## add specified number of trailing spaces to variable names:
          plotnames <- varimps$varname
          if (is.numeric(abbreviate) && abbreviate > 0) {
            plotnames <- abbreviate(plotnames, minlength = abbreviate)
          }
          if (diag.xlab.vert > 0) {
            for (i in 1:diag.xlab.vert) {
              plotnames <- paste0(plotnames, " ")
            }
          }
          text(xlab.pos + diag.xlab.hor, par("usr")[3], srt = 45, adj = 1, xpd = TRUE, 
               labels = plotnames, cex = cex.axis)
        } else {
          plotnames <- varimps$varname
          if (is.numeric(abbreviate) && abbreviate > 0) {
            plotnames <- abbreviate(plotnames, minlength = abbreviate)
          }
          barplot(height = varimps$imp, names.arg = plotnames, ylab = ylab,
                main = main, ...)
        }
      }
    }
    
    if (!is.na(round)) {
      baseimps[,sapply(baseimps, is.numeric)] <- round(baseimps[,sapply(baseimps, is.numeric)], digits = round)
      varimps[,sapply(varimps, is.numeric)] <- round(varimps[,sapply(varimps, is.numeric)], digits = round)
    }
    
    if (object$family %in% c("mgaussian","multinomial")) {
      keep <- c("rule", "description", gsub("coefficient", "importance", coef_inds), 
                coef_inds, "sd")
    } else {
      keep <- c("rule", "description", "imp", "coefficient", "sd")
    }
    
    baseimps <- data.frame(baseimps[, keep], stringsAsFactors = FALSE)
    row.names(baseimps) <- row.names(varimps) <- NULL
    
    return(invisible(list(varimps = varimps, baseimps = baseimps)))
    
  } else {
    warning("No non-zero terms in the ensemble. All importances are zero.")
    return(invisible(NULL))
  }
}





#' Compute bootstrapped null interaction prediction rule ensembles
#'
#' \code{bsnullinteract} generates bootstrapped null interaction models,
#' which can be used to derive a reference distribution of the test statistic
#' calculated with \code{\link{interact}}.
#'
#' @param object object of class \code{\link{pre}}.
#' @param nsamp numeric. Number of bootstrapped null interaction models to be
#' derived.
#' @inheritParams print.pre
#' @param parallel logical. Should parallel foreach be used to generate initial
#' ensemble? Must register parallel beforehand, such as doMC or others.
#' @param verbose logical. should progress be printed to the command line?
#' @return A list of length \code{nsamp} with null interaction models, to be
#' used as input for \code{\link{interact}}.
#' @examples \donttest{set.seed(42)
#' airq.ens <- pre(Ozone ~ ., data=airquality[complete.cases(airquality),])
#' nullmods <- bsnullinteract(airq.ens)
#' interact(airq.ens, nullmods = nullmods, col = c("#7FBFF5", "#8CC876"))}
#' @details Note that computation of bootstrapped null interaction models is 
#' computationally intensive. The default number of samples is set to 10,
#' but for reliable results argument \code{nsamp} should be set to a higher
#' value (e.g., \eqn{\ge 100}). 
#' 
#' See also section 8.3 of Friedman & Popescu (2008).
#' 
#' @references Fokkema, M. (2018). Fitting prediction rule ensembles with R package pre.
#' \url{https://arxiv.org/abs/1707.07149}.
#' 
#' Friedman, J. H., & Popescu, B. E. (2008). Predictive learning 
#' via rule ensembles. \emph{The Annals of Applied Statistics, 2}(3), 916-954.
#' @seealso \code{\link{pre}}, \code{\link{interact}} 
#' @export
bsnullinteract <- function(object, nsamp = 10, parallel = FALSE,
                           penalty.par.val = "lambda.1se", verbose = FALSE)
{
  
  ## TODO: Implement for all types of response variables
  if (object$family %in% c("mgaussian", "multinomial", "cox", "binomial")) {
    stop("Function bsnullinteract not implemented yet for binomial, multinomial and multivariate outcomes.")
  }
  
  if(parallel) {
    if (!(requireNamespace("foreach"))) {
      warning("Parallel computation of function bsnullinteract() requires package foreach,
              which is currently not installed. Argument parallel will be set to FALSE.
              To run in parallel, download and install package foreach from CRAN, and run again.")
      parallel <- FALSE
    }
  }
  # create call for generating bootstrapped null models:
  bsnullmodcall <- object$call
  bsnullmodcall$maxdepth <- 1L
  # create call for model allowing for interactions, grown on bootstrapped
  # datasets without interactions:
  bsintmodcall <- object$call
  bsintmodcall$verbose <- FALSE
  # compute bootstrapped null datasets (i.e., datasets with no interactions):
  if (parallel) {
    if (verbose) cat("This may take a while.")
    bs.ens <- foreach::foreach(i = 1:nsamp, .packages = "pre") %dopar% {
      # step 1: Take bootstrap sample {x_p, y_p}:
      bs_inds <- sample(1:nrow(object$data), nrow(object$data), replace = TRUE)
      bsdataset <- object$data[bs_inds,]
      # step 2: Build F_A, a null interaction model involving main effects only using {x_p, y_p}:
      bsnullmodcall$data <- bsdataset
      bs.ens.null <- eval(bsnullmodcall)
      # step 3: first part of formula 47 of F&P2008:
      # Calculate predictions F_A(x) for original x, using the null interaction model F_A:
      F_A_of_x <- predict.pre(bs.ens.null, newdata = object$data, 
                              penalty.par.val = penalty.par.val)
      # step 4: third part of formula 47 of F&P2008:
      # Calculate predictions F_A(x_p):
      F_A_of_x_p <- predict.pre(bs.ens.null, newdata = bsdataset,
                                penalty.par.val = penalty.par.val)
      # step 5: Calculate ytilde of formula 47 of F&P2008:
      ytilde <- F_A_of_x + object$data[bs_inds, object$y_names] - F_A_of_x_p
      # step 6: Build a model using (x,ytilde), using the same procedure as was
      # originally applied to (x,y):
      bsintmodcall$data <- object$data
      bsintmodcall$data[,all.vars(object$call$formula[[2]])] <- ytilde
      eval(bsintmodcall)
    }
  } else {
    bs.ens <- list()
    if (verbose) cat("This may take a while. Computing null model ")
    for(i in 1:nsamp) {
      if (verbose) {cat(i, "of", nsamp, ", ")}
      # step 1: Take bootstrap sample {x_p, y_p}:
      bs_inds <- sample(1:nrow(object$data), nrow(object$data), replace = TRUE)
      bsdataset <- object$data[bs_inds,]
      # step 2: Build F_A, a null interaction model involving main effects only using {x_p, y_p}:
      bsnullmodcall$data <- as.symbol(quote(bsdataset))
      bs.ens.null <- eval(bsnullmodcall, envir = environment())
      # step 3: first part of formula 47 of F&P2008:
      # Calculate predictions F_A(x) for original x, using the null interaction model F_A:
      F_A_of_x <- predict.pre(bs.ens.null, newdata = object$data, type = "link")
      # step 4: third part of formula 47 of F&P2008:
      # Calculate predictions F_A(x_p):
      F_A_of_x_p <- predict.pre(bs.ens.null, newdata = bsdataset,
                                penalty.par.val = penalty.par.val,
                                type = "link")
      # step 5: Calculate ytilde of formula 47 of F&P2008:
      # TODO: Does not compute for categorical outcomes: 
      ytilde <- F_A_of_x + (object$data[bs_inds, object$y_names] - F_A_of_x_p)
      # step 6: Build a model using (x,ytilde), using the same procedure as was
      # originally applied to (x,y):
      tmp <- object$data
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
  preds_x <- predict.pre(object, newdata = object$data, penalty.par.val = penalty.par.val)
  # Calculate the expected value of F_j(x_j), over all observed values x_/j,
  # and the expected value of F_/j(x_/j), over all observed values x_j:
  exp_dataset <- object$data[rep(row.names(object$data),
                                      times = nrow(object$data)),]
  exp_dataset[,varname] <- rep(object$data[,varname], each = nrow(object$data))
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
  exp_dataset$i_xj <- rep(1:nrow(object$data), each = nrow(object$data))
  preds_xj <- aggregate(yhat ~ i_xj, data = exp_dataset, FUN = mean)$yhat
  # expected value of F_/j(x_/j), over all observed values x_j:
  exp_dataset$i_xnotj <-  rep(1:nrow(object$data), times = nrow(object$data))
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





#' Calculate interaction statistics for variables in a prediction rule ensemble 
#' (pre)
#'
#' \code{interact} calculates test statistics for assessing the strength of
#' interactions between a set of user-specified input variable(s), and all 
#' other input variables.
#'
#' @param object an object of class \code{\link{pre}}.
#' @param varnames character vector. Names of variables for which interaction
#' statistics should be calculated. If \code{NULL}, interaction statistics for
#' all predictor variables with non-zeor coefficients will be calculated (which
#' may take a long time).
#' @param nullmods object with bootstrapped null interaction models, resulting
#' from application of \code{bsnullinteract}.
#' @inheritParams print.pre
#' @param quantprobs numeric vector of length two. Probabilities that should be
#' used for plotting the range of bootstrapped null interaction model statistics.
#' Only used when \code{nullmods} argument is specified and \code{plot = TRUE}.
#' The default yields sample quantiles corresponding to .05 and .95 probabilities.  
#' @param plot logical. Should interaction statistics be plotted?
#' @param col character vector of length one or two. The first value specifies 
#' the color to be used for plotting the interaction statistic from the training
#' data, the second color is used for plotting the interaction statistic from 
#' the bootstrapped null interaction models. Only used when \code{plot = TRUE}. 
#' Only the first element will be used if \code{nullmods = NULL}.
#' @param ylab character string. Label to be used for plotting y-axis.
#' @param legend.text character vector of length two to be used for plotting
#' the legend. Only used when \code{nullmods} is specified. If \code{FALSE},
#' no legend is plotted.
#' @param main character. Main title for the bar plot.
#' @param  se.linewidth numeric. Width of the whiskers of the plotted standard 
#' error bars (in inches).
#' @param k integer. Calculating interaction test statistics is computationally
#' intensive, so  calculations are split up in several parts to prevent memory
#' allocation errors. If a memory allocation error still occurs, increase k.
#' @param verbose logical. Should progress information be printed to the
#' command line?
#' @param parallel logical. Should parallel foreach be used? Must register
#' parallel beforehand, such as doMC or others.
#' @param ... Additional arguments to be passed to \code{barplot}.
#' @examples
#' \donttest{set.seed(42)
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
#' have not yet been studied in detail, but results are likely to get more 
#' reliable when the number of bootstrapped null interaction models is larger.
#' The default of the \code{bsnullinteract} function is to generate 10 
#' bootstrapped null interaction datasets, to yield shorter computation times.
#' To obtain a more reliable result, however, users are advised to
#' set the \code{nsamp} argument \eqn{\ge 100}.
#' 
#' See also section 8 of Friedman & Popescu (2008).
#' 
#' @references Fokkema, M. (2018). Fitting prediction rule ensembles with R package pre.
#' \url{https://arxiv.org/abs/1707.07149}.
#' 
#' Friedman, J. H., & Popescu, B. E. (2008). Predictive learning 
#' via rule ensembles. \emph{The Annals of Applied Statistics, 2}(3), 916-954.
#' @seealso \code{\link{pre}}, \code{\link{bsnullinteract}} 
#' @export
interact <- function(object, varnames = NULL, nullmods = NULL, 
                     penalty.par.val = "lambda.1se", quantprobs = c(.05, .95),
                     plot = TRUE, col = c("darkgrey", "lightgrey"), 
                     ylab = "Interaction strength", 
                     main = "Interaction test statistics", 
                     se.linewidth = .05, 
                     legend.text = c("observed", "null model median"),
                     parallel = FALSE, k = 10, verbose = FALSE, ...) {
  
  ## TODO: Implement for all response variable types
  if (object$family %in% c("mgaussian", "multinomial")) {
    stop("Function interact not implemented yet for multivariate and multinomial outcomes.")
  }

  # Preliminaries:
  if(parallel) {
    if (!(requireNamespace("foreach"))) {
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
    H <- foreach::foreach(i = 1:length(varnames), .combine = "c", .packages = "pre") %dopar% {
      # Calculate H_j for the original dataset:
      Hsquaredj(object = object, varname = varnames[i], k = k,
                penalty.par.val = penalty.par.val, verbose = verbose)
    }
    names(H) <- varnames
    if (!is.null(nullmods)) {
      nullH <- foreach::foreach(i = 1:length(varnames), .combine = "cbind", .packages = "pre") %dopar% {
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
      observed <- rbind(H, apply(nullH, 2, mean))
      H0_medians <- apply(nullH, 2, median)
      lower_quant <- apply(nullH, 2, quantile, probs = quantprobs[1])
      upper_quant <- apply(nullH, 2, quantile, probs = quantprobs[2])
      x_coords <- barplot(observed, beside = TRUE, 
                          ylim = c(0, max(upper_quant, observed)), 
                          las = 1, main = main, col = col, 
                          legend.text = legend.text, ...)
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
#' @inheritParams print.pre
#' @param linear.terms logical. Should linear terms be included in the plot?
#' @param nterms numeric. The total number of terms (or rules, if 
#' \code{linear.terms = FALSE}) being plotted. Default is \code{NULL}, 
#' resulting in all terms of the final ensemble to be plotted.
#' @param fill character of length 1 or 2. Background color(s) for terminal 
#' panels. If one color is specified, all terminal panels will have the 
#' specified background color. If two colors are specified (the default, the 
#' first color will be used as the background color for rules with a positively 
#' valued coefficient; the second color for rules with a negatively valued 
#' coefficient.
#' @param plot.dim integer vector of length two. Specifies the number of rows
#' and columns in the plot. The default yields a plot with three rows and three 
#' columns, depicting nine baselearners per plotting page.
#' @param ask logical. Should user be prompted before starting a new page of
#' plots?
#' @param exit.label character string. Label to be printed in nodes to which 
#' the rule does not apply (``exit nodes'')?
#' @param standardize logical. Should printed importances be standardized? See
#' \code{\link{importance}}.
#' @param ... Arguments to be passed to \code{\link[grid]{gpar}}.
#' @examples
#' \donttest{set.seed(42)
#'  airq.ens <- pre(Ozone ~ ., data = airquality[complete.cases(airquality),])
#'  plot(airq.ens)}
#' @seealso \code{\link{pre}}, \code{\link{print.pre}}
#' @method plot pre
#' @export
plot.pre <- function(x, penalty.par.val = "lambda.1se", linear.terms = TRUE, 
                     nterms = NULL, fill = "white", ask = FALSE, 
                     exit.label = "0", standardize = FALSE, plot.dim = c(3, 3), 
                     ...) {
  
  ## rpart uses < and >=, whereas partykit uses <= and > for splits. 
  ## This should be supplied to partysplit for plotting:
  if (is.null(x$call$tree.unbiased)) {
    right <- TRUE
  } else if (x$call$tree.unbiased) {
    right <- TRUE  
  } else if (!x$call$tree.unbiased) {
    right <- FALSE
  }
  
  if (x$family %in% c("mgaussian", "multinomial")) {
    warning("Plotting function not yet fully functional for multivariate and multinomial outcomes.")
  }
  
  if (!(requireNamespace("grid"))) {
    stop("Function plot.pre requires package grid. Download and install package grid from CRAN, and run again.")
  }

  ## Get nonzero terms:
  if (x$family %in% c("multinomial", "mgaussian")) {
    coefs <- coef(x, penalty.par.val = penalty.par.val)
    nonzeroterms <- coefs[rowSums(coefs[,!names(coefs) %in% c("rule", "description")]) != 0,]
    if ("(Intercept)" %in% nonzeroterms$rule) {
      intercept <- nonzeroterms[which(nonzeroterms$rule == "(Intercept)"), "coefficient"] # may be needed for plotting linear terms later      
      nonzeroterms <- nonzeroterms[-which(nonzeroterms$rule == "(Intercept)"), ] # omit intercept
    }
  } else {
    nonzeroterms <- importance(x, plot = FALSE, global = TRUE, 
                               penalty.par.val = penalty.par.val, 
                               standardize = standardize)$baseimps
  }

  if (!linear.terms) {
    nonzeroterms <- nonzeroterms[grep("rule", nonzeroterms$rule),]
  }
  if (!is.null(nterms) && nrow(nonzeroterms) > nterms) {
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
  n_terms_per_plot <- plot.dim[1L] * plot.dim[2L]
  nplots <- ceiling(nrow(nonzeroterms) / n_terms_per_plot)
  nonzeroterms$plotno <- rep(1L:nplots, each = n_terms_per_plot)[1L:nrow(nonzeroterms)]
  nonzeroterms$rowno <- rep(rep(1L:plot.dim[1L], each = plot.dim[2L]), length.out = nrow(nonzeroterms))
  nonzeroterms$colno <- rep(rep(1L:plot.dim[2L], times = plot.dim[1L]), length.out = nrow(nonzeroterms))
  
  ## Generate a plot for every baselearner:
  for(i in 1:nrow(nonzeroterms)) {
    
    if (conditions[[i]][1L] == "linear") { 
      ## Plot linear term:
      ## Open new plotting page if needed:
      if (nonzeroterms$rowno[i] == 1L && nonzeroterms$colno[i] == 1L) {
        grid::grid.newpage()
        grid::pushViewport(grid::viewport(layout = grid::grid.layout(plot.dim[1L], plot.dim[2L])))
      }
      ## open correct viewport:
      grid::pushViewport(grid::viewport(layout.pos.col = nonzeroterms$colno[i],
                                        layout.pos.row = nonzeroterms$rowno[i]))
      ## Plot the linear term:

      if (x$family %in% c("mgaussian", "multinomial")) {
        coef_names <- names(nonzeroterms)[grepl("coefficient.", names(nonzeroterms))]        
        coef_names <- data.frame(name = coef_names,
                         value = t(round(nonzeroterms[i,coef_names], digits = 3L)))
        coef_names <- paste0(apply(coef_names, 1L, paste0, collapse = " = "), collapse = "\n")
        grid::grid.text(paste0("Linear effect of ", nonzeroterms$rule[i], 
                               "\n\n", coef_names), gp = grid::gpar(...))
      } else {
        ## This seems to work for plotting but should be tested::
        #lattice::xyplot(y ~ x, 
        #                data = data.frame(y = range(x$data[,x$y_names]), x = range(x$data[,nonzeroterms[i, "rule"]])), 
        #                type = "n", ylab = x$y_names, xlab = nonzeroterms[i, "rule"], main = paste("Linear effect of", nonzeroterms$rule[i]),
        #                panel = function(...) {
        #                  lattice::panel.abline(a = intercept, b = nonzeroterms[i, "coefficient"])
        #                  lattice::panel.xyplot(...)
        #                })
        grid::grid.text(paste0("Linear effect of ", nonzeroterms$rule[i], 
                             "\n\n Coefficient = ", round(nonzeroterms$coefficient[i], digits = 3L),
                               "\n\n Importance = ", round(nonzeroterms$imp[i], digits = 3L)),
                        gp = grid::gpar(...))        
      }
      grid::popViewport()
      
    } else { ## Otherwise, plot rule:
      
      # Create lists of arguments and operators for every condition:
      cond <- list()
      # check whether the operator is " < ", " <= " or "%in%  "
      # split the string using the operator, into the variable name and splitting value, 
      # which is used to define split = partysplit(id, value)
      # make it a list:
      for (j in 1L:length(conditions[[i]])) {
        ## TODO: see get_conditions() function below for possible improvements to this code:
        condition_j <- conditions[[i]][[j]]
        cond[[j]] <- character()
        if (length(grep(" > ", condition_j)) > 0) {
          cond[[j]][1L] <- unlist(strsplit(condition_j, " > "))[1L]
          cond[[j]][2L] <- " > "
          cond[[j]][3L] <- unlist(strsplit(condition_j, " > "))[2L]
        } else if (length(grep(" >= ", condition_j)) > 0) {
          cond[[j]][1L] <- unlist(strsplit(condition_j, " >= "))[1L]
          cond[[j]][2L] <- " >= "
          cond[[j]][3L] <- unlist(strsplit(condition_j, " >= "))[2L]
        } else if (length(grep(" <= ", condition_j)) > 0) {
          cond[[j]][1L] <- unlist(strsplit(condition_j, " <= "))[1L]
          cond[[j]][2L] <- " <= "
          cond[[j]][3L] <- unlist(strsplit(condition_j, " <= "))[2L]
        } else if (length(grep(" < ", condition_j)) > 0) {
          cond[[j]][1L] <- unlist(strsplit(condition_j, " < "))[1L]
          cond[[j]][2L] <- " < "
          cond[[j]][3L] <- unlist(strsplit(condition_j, " < "))[2L]
        } else if (length(grep(" %in% ", condition_j)) > 0) {
          cond[[j]][1L] <- unlist(strsplit(condition_j, " %in% "))[1L]
          cond[[j]][2L] <- " %in% "
          cond[[j]][3L] <- unlist(strsplit(condition_j, " %in% "))[2L]
        }
      }
      ncond <- length(cond)
      cond <- rev(cond)
      
      ## generate empty dataset for all the variables appearing in the rules: 
      treeplotdata <- data.frame(matrix(ncol = ncond))
      for (j in 1L:ncond) {
        names(treeplotdata)[j] <- cond[[j]][1L]
        if (cond[[j]][2L] == " %in% ") {
          treeplotdata[ , j] <- factor(treeplotdata[ , j])
          faclevels <- substring(cond[[j]][3L], first = 2L)
          faclevels <- gsub(pattern = "\"", replacement = "", x = faclevels, fixed = TRUE)
          faclevels <- gsub(pattern = "(", replacement = "", x = faclevels, fixed = TRUE)
          faclevels <- gsub(pattern = ")", replacement = "", x = faclevels, fixed = TRUE)
          faclevels <- unlist(strsplit(faclevels, ", ",))
          levels(treeplotdata[,j]) <- c(
            levels(x$data[ , cond[[j]][1L]])[levels(x$data[ , cond[[j]][1L]]) %in% faclevels],
            levels(x$data[ , cond[[j]][1L]])[!(levels(x$data[ , cond[[j]][1L]]) %in% faclevels)])
          cond[[j]][3L] <- length(faclevels)
        }
      }

      
      ## Generate partynode objects for plotting:
      
      nodes <- list()
      ## Create level 0 (bottom level, last two nodes):
      ## If last condition has " > " : exit node on left, coefficient on right:
      if (cond[[1L]][2L] %in% c(" > ", " >= ")) { # If condition involves " > ", the tree has nonzero coef on right:
        nodes[[1L]] <- list(id = 1L, split = NULL, kids = NULL, surrogates = NULL, 
                           info = exit.label)
        if (x$family %in% c("multinomial", "mgaussian")) {
          info <- paste(round(nonzeroterms[i, grep("coefficient", names(nonzeroterms))], digits = 3L), collapse = "\n")
          nodes[[2L]] <- list(id = 2L, split = NULL, kids = NULL, surrogates = NULL,
                             info = info)
        } else {
          nodes[[2L]] <- list(id = 2L, split = NULL, kids = NULL, surrogates = NULL,
                             info = round(nonzeroterms$coefficient[i], digits = 3L))  
        }
      } else { 
        ## If last condition has " <= " or " %in% " : coefficient on left, exit node on right:
        if (x$family %in% c("multinomial", "mgaussian")) {
          info <- paste(round(nonzeroterms[i, grep("coefficient", names(nonzeroterms))], digits = 3L), collapse = "\n")
          nodes[[1L]] <- list(id = 1L, split = NULL, kids = NULL, surrogates = NULL,
                             info = info)
        } else {
          nodes[[1L]] <- list(id = 1L, split = NULL, kids = NULL, surrogates = NULL,
                             info = round(nonzeroterms$coefficient[i], digits = 3L))
        }
        nodes[[2]] <- list(id = 2L, split = NULL, kids = NULL, surrogates = NULL,
                           info = exit.label)
      }
      class(nodes[[1L]]) <- class(nodes[[2L]]) <- "partynode"
      
      ## Create inner levels (if necessary):
      if (ncond > 1L) {
        for (level in 1L:(ncond - 1L)) {
          if (cond[[level + 1L]][2L] == " > ") { 
            ## If condition in level above has " > " : exit node on left, right node has kids:
            nodes[[level * 2L + 1L]] <- list(id = as.integer(level * 2L + 1L), 
                                          split = NULL, 
                                          kids = NULL, 
                                          surrogates = NULL, 
                                          info = exit.label)
            nodes[[level * 2L + 2L]] <- list(id = as.integer(level * 2L + 2L), 
                                          split = partysplit(as.integer(level), 
                                                             breaks = as.numeric(cond[[level]][3L]),
                                                             right = right),
                                          kids = list(nodes[[level * 2L - 1L]], nodes[[level * 2L]]),
                                          surrogates = NULL, 
                                          info = NULL)
          } else { 
            ## If condition in level above has " <= " or " %in% " : left node has kids, exit node right:
            nodes[[level * 2L + 1L]] <- list(id = as.integer(level * 2L + 1L),
                                          split = partysplit(as.integer(level), 
                                                             breaks = as.numeric(cond[[level]][3]),
                                                             right = right),
                                          kids = list(nodes[[level * 2 - 1]], nodes[[level * 2]]),
                                          surrogates = NULL, 
                                          info = NULL)
            nodes[[level * 2L + 2L]] <- list(id = as.integer(level * 2L + 2L), 
                                          split = NULL,
                                          kids = NULL, 
                                          surrogates = NULL, 
                                          info = exit.label)
          }  
          class(nodes[[level * 2L + 1L]]) <- class(nodes[[level * 2L + 2L]]) <- "partynode"
        }
      }
      
      ## Create root node:
      nodes[[ncond * 2L + 1L]] <- list(id = as.integer(ncond * 2L + 1L),
                                     split = partysplit(as.integer(ncond), 
                                                        breaks = as.numeric(cond[[ncond]][3L]),
                                                        right = right),
                                     kids = list(nodes[[ncond * 2L - 1L]], nodes[[ncond * 2L]]),
                                     surrogates = NULL, 
                                     info = NULL)
      class(nodes[[ncond * 2L + 1L]]) <- "partynode"
      
    
      ## Open new plotting page if needed:
      if (nonzeroterms$rowno[i] == 1L && nonzeroterms$colno[i] == 1L) {
        grid::grid.newpage()
        grid::pushViewport(grid::viewport(layout = grid::grid.layout(plot.dim[1L], plot.dim[2L])))
      }
    
    
      ## Open viewport:
      grid::pushViewport(grid::viewport(layout.pos.col = nonzeroterms$colno[i],
                                        layout.pos.row = nonzeroterms$rowno[i]))
    
      ## Plot the rule:
      fftree <- party(nodes[[ncond * 2L + 1L]], data = treeplotdata)
      if (x$family %in% c("mgaussian", "multinomial")) {
        if (x$family == "mgaussian") {
          ht <- length(x$y_names)
        } else {
          ht <- nlevels(x$data[,x$y_names])
        }
        plot(fftree, newpage = FALSE, main = nonzeroterms$rule[i],
             inner_panel = node_inner(fftree, id = FALSE),
             terminal_panel = node_terminal(fftree, id = FALSE, 
                                            fill = "white", height = ht), 
             gp = grid::gpar(...))
      } else {
        plot(fftree, newpage = FALSE, 
             main = paste0(nonzeroterms$rule[i], ": Importance = ", round(nonzeroterms$imp[i], digits = 3L)),
             inner_panel = node_inner(fftree, id = FALSE),
             terminal_panel = node_terminal(fftree, id = FALSE, fill = ifelse(length(fill) > 1L, 
                                                                              ifelse(nonzeroterms$coefficient[i] > 0, fill[1L], fill[2L]), 
                                                                              fill)), 
             gp = grid::gpar(...))      
      }
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



## Get rule conditions
## 
## \code{get_conditions} returns rule conditions in a matrix form
##  
## @param object object of class pre
## @param penalty.par.val character. Value of the penalty parameter value 
## \eqn{\lambda} to be used for selecting the final ensemble. The ensemble 
## with penalty parameter criterion yielding minimum cv error 
## (\code{"lambda.min"}) is taken, by default. Alternatively, the penalty 
## parameter yielding error within 1 standard error of minimum cv error 
## ("\code{lambda.1se}"), or a numeric value may be specified, corresponding 
## to one of the values of lambda in the sequence used by glmnet,
## for which estimated cv error can be inspected by running \code{x$glmnet.fit}
## and \code{plot(x$glmnet.fit)}.
## @examples \donttest{set.seed(42)
## airq.ens <- pre(Ozone ~ ., data = airquality)
## get_conditions(airq.ens)}
get_conditions <- function(object, penalty.par.val = "lambda.1se") {
  ## get maximum rule depth used for generating ensemble:
  if (is.null(object$call$maxdepth)) {
    maxdepth <- 3
  } else {
    maxdepth <- object$call$maxdepth
  }
  ## get the rules from the ensemble:
  rules <- object$rules  
  ## cut rules into parts:
  parts <- strsplit(rules[,"description"], split = " ")
  ## turn parts into matrix: 
  parts <- matrix(unlist(lapply(parts, `length<-`, max(lengths(parts)))), 
                  ncol = max(lengths(parts)), byrow = TRUE)
  ## eliminate every fourth column (has "&" only):
  parts <- data.frame(parts[,!(1:ncol(parts) %% 4 == 0)])
  ## set every third column to type numeric (not a good idea, can be factors, too):
  parts[,(1:ncol(parts) %% 3 == 0)] <- apply(parts[,(1:ncol(parts) %% 3 == 0)], 2, as.numeric)
  parts <- data.frame(rule = rules$rule, parts)
  names(parts) <- c("rule", paste0(rep(c("splitvar", "splitop", "splitval"), 
                                       times = maxdepth), 
                                   rep(1:maxdepth, each = 3)))
  ## only get rules that are in final ensemble:
  coefs <- coef(object, penalty.par.val = penalty.par.val)
  nonzero_rules <- coefs[coefs$coefficient != 0,]$rule
  parts <- parts[parts$rule %in% nonzero_rules,]
  ## return result:
  return(parts)
}






#' Plot correlations between baselearners in a prediction rule ensemble (pre)
#' 
#' \code{corplot} plots correlations between baselearners in a prediction rule ensemble
#'  
#' @param object object of class pre
#' @inheritParams print.pre
#' @param colors vector of contiguous colors to be used for plotting. If 
#' \code{colors = NULL} (default), \code{colorRampPalette} is used to generate
#' a sequence of 200 colors going from red to white to blue. A different set of 
#' plotting colors can be specified here, for example: 
#' \code{cm.colors(100)}, \code{colorspace::rainbow_hcl)(100)} 
#' or \code{colorRampPalette(c("red", "yellow", "green"))(100)}.
#' @param fig.plot plotting region to be used for correlation plot. See 
#' \code{fig} under \code{\link{par}}.
#' @param fig.legend plotting region to be used for legend. See \code{fig} 
#' under \code{\link{par}}.
#' @param legend.breaks numeric vector of breakpoints to be depicted in the 
#' plot's legend. Should be a sequence from -1 to 1.
#' @examples \donttest{set.seed(42)
#' airq.ens <- pre(Ozone ~ ., data = airquality[complete.cases(airquality),])
#' corplot(airq.ens)}
#' @seealso See
#' \code{\link[colorspace]{rainbow_hcl}} and \code{\link[grDevices]{colorRampPalette}}.
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
        col = colors)
  axis(1, at = 1:nrow(cormat), labels = colnames(cormat), las = 2)
  axis(2, at = 1:ncol(cormat), labels = colnames(cormat), las = 2)
  ## plot legend:
  par(fig = fig.legend, new = TRUE)
  col_inds <- round(seq(1, length(colors), length.out = length(legend.breaks)-1))
  image.scale(z = legend.breaks,
              col = colors[col_inds],
              axis.pos = 4,
              breaks = legend.breaks)
  axis(4, at = legend.breaks, las = 2)
  return(invisible(cormat))
}



#' Explain predictions from final prediction rule ensemble
#'
#' \code{explain} shows which rules apply to which observations and visualizes
#' the contribution of rules and linear predictors to the predicted values
#'
#' @param object object of class \code{\link{pre}}.
#' @param newdata optional dataframe of new (test) observations, including all
#' predictor variables used for deriving the prediction rule ensemble.
#' @inheritParams print.pre
#' @param plot logical. Should explanations be plotted?
#' @param intercept logical. Specifies whether intercept should be included in
#' explaining predictions.
#' @param center.linear logical. Specifies whether linear terms should be
#' centered with respect to the training sample mean before computing their 
#' contribution to the predicted value.
#' @param pred.type character. Specifies the type of predicted values provided
#' in the plot(s).
#' @param plot.max.nobs numeric. Specifies maximum number of observations
#' for which explanations will be plotted. The default \code{4} plots the
#' explanation for the first four observations supplied in \code{newdata}.
#' @param plot.dim numeric vector of length 2. Specifies the number of rows and
#' columns in the resulting plot.
#' @param cex numeric. Specifies the relative text size of title, tick and axis
#' labels.
#' @param plot.obs.names logical vector of length 1, NULL, or character vector 
#' of length \code{nrow(data)} supplying the names that should be used for
#' individual observations' plots. If \code{TRUE} (default), 
#' \code{rownames(newdata)} will be used as titles. If \code{NULL},
#' \code{paste("Observation", 1:nrow(newdata))} will be used as titles. If 
#' \code{FALSE}, no titles will be plotted.
#' @param digits integer. Specifies the number of digits used in depcting the
#' predicted values in the plot.
#' @param bar.col character vector of length two. Specifies the colors to be used for
#' plotting the positive and negative contributions to the predictions, respectively.
#' @param rule.col character. Specifies the color to be used for plotting the rule
#' descriptions. If \code{NULL}, rule descriptions are not plotted.
#' @param ylab character. Specifies the label for the horizonantal (y-) axis.
#' @details Provides a graphical depiction of the contribution of rules and
#' linear terms to the individual predictions (if \code{plot = TRUE}.
#' Invisibly returns a list with objects \code{predictors} and
#' \code{contribution}. \code{predictors} contains the values of the rules and
#' linear terms for each observation in \code{newdata}, for those rules
#' and linear terms included in the final ensemble with the specified
#' value of \code{penalty.par.val}. \code{contribution} contains the
#' values of \code{predictors}, multiplied by the estimated values
#' of the coefficients in the final ensemble selected with the
#' specified value of \code{penalty.par.val}.
#' All contributions are calculated w.r.t. the intercept, by default.
#' Thus, if a given rule applies to an observation in \code{newdata}, 
#' the contribution of that rule equals the estimated coefficient of
#' that rule. If a given rule does not apply to an observation in
#' \code{newdata}, the contribution of that rule equals 0. 
#' For linear terms, contributions can be centered, or not (the default).
#' Thus, by default the contribution of a linear terms for an 
#' observation in \code{newdata} equals the obeservation's value of the 
#' linear term, times the estimated coefficient of the linear term.
#' If \code{center.linear = TRUE}, the contribution of a linear term 
#' for an observation in \code{newdata} equals (the value of the linear
#' temr, minus the mean value of the linear term in the training data)
#' times the estimated coefficient for the linear term.
#' @examples \donttest{airq <- airquality[complete.cases(airquality), ]
#' set.seed(1)
#' train <- sample(1:nrow(airq), size = 100)
#' set.seed(42)
#' airq.ens <- pre(Ozone ~ ., data = airq[train,])
#' airq.ens.exp <- explain(airq.ens, newdata = airq[-train,])
#' airq.ens.exp$predictors
#' airq.ens.exp$contribution
#' 
#' ## Fit PRE with three trees and lower penalty parameter value,
#' ## to allow linear terms to appear in final ensemble:
#' set.seed(42)
#' airq.ens2 <- pre(Ozone ~ ., data = airq[train,], ntrees = 3L)
#' ## Compared to intercept, Month has negative and 
#' ## Day has positive contribution:
#' explain(airq.ens2, newdata = airq[-train,][1:2,], 
#'         penalty.par.val = "lambda.min")$contribution
#' ## Compared with training data means, Month has positive 
#' ## and Day has negative contribution: 
#' explain(airq.ens2, newdata = airq[-train,][1:2,], 
#'         penalty.par.val = "lambda.min", center.linear = TRUE)$contribution
#' }
#' @seealso \code{\link{pre}}, \code{\link{plot.pre}},
#' \code{\link{coef.pre}}, \code{\link{importance}}, \code{\link{cvpre}},
#' \code{\link{interact}}, \code{\link{print.pre}}
#' @export
explain <- function(object, newdata, penalty.par.val = "lambda.1se",
                    plot = TRUE, intercept = FALSE, center.linear = FALSE,
                    plot.max.nobs = 4, plot.dim = c(2, 2),
                    plot.obs.names = TRUE, pred.type = "response", 
                    digits = 3L,  cex = .8,
                    ylab = "Contribution to linear predictor",
                    bar.col = c("#E495A5", "#39BEB1"),
                    rule.col = "darkgrey") {
  
  ## check arguments:
  if (!inherits(object, what = "pre")) {
    stop("Argument object should specfify an object of class 'pre'.")
  }
  if (!is.data.frame(newdata)) {
    stop("Argument newdata should specify a data frame.")
  }
  if (!(is.logical(intercept) || length(intercept) != 1L)) {
    stop("Argument intercept should be a logical vector of length 1.")
  }
  if (!(is.logical(plot) || length(plot) != 1L)) {
    stop("Argument plot should be a logical vector of length 1.")
  }
  
  preds <- round(predict(object, newdata = newdata, type = pred.type,
                         penalty.par.val = penalty.par.val), digits = digits)
  
  ## Prepare newdata:
  winsfrac <- (object$call)$winsfrac
  if (is.null(winsfrac)) winsfrac <- formals(pre)$winsfrac
  
  ## Check if variable names and classes are the same in newdata as in object$data:
  if (!all(object$x_names %in% names(newdata))) {
    newdata <- model.frame(as.Formula((object$call)$formula), data = newdata, 
                           rhs = NULL, lhs = 0, na.action = NULL)
  } else {
    newdata <- newdata[ , object$x_names]
  }
  
  ## Coerce character and logical variables to factors:
  if (any(char_names <- sapply(newdata, is.character))) {
    char_names <- names(newdata)[char_names]
    data[ , char_names] <- sapply(newdata[ , char_names], factor)
  }
  if (any(logic_names <- sapply(newdata, is.logical))) {
    logic_names <- names(newdata)[logic_names]
    newdata[ , logic_names] <- sapply(newdata[ , logic_names], factor)
  } 
  
  ## Coerce ordered categorical variables to numeric, if necessary:
  if (if (is.null((object$call)$ordinal)) {
    formals(pre)$ordinal
  } else {
    (object$call)$ordinal
  }) {
    if (any(ord_var_inds <- sapply(newdata, is.ordered))) {
      newdata[ , ord_var_inds] <- sapply(newdata[ , ord_var_inds], as.numeric)
    }
  }
  
  if (any(is.na(newdata))) {
    newdata <- newdata[complete.cases(newdata),]
    warning("Some observations in newdata have missing predictor variable values and will be removed.", immediate. = TRUE)
  }
  
  ## Check and set factor levels of newdata to variable levels in object$data:
  if (any(factor_inds <- sapply(newdata, is.factor))) {
    for (i in names(newdata)[factor_inds]) {
      if (all(levels(newdata[ , i]) %in% levels(object$data[ , i]))) {
        levels(newdata[ , i]) <- levels(object$data[ , i])
      } else {
        stop("Variable ", i, " has levels not present in training data. Cannot compute predictions.")
      }
    }
  }
  
  newdata <- get_modmat(
    wins_points = object$wins_points, 
    x_scales = object$x_scales, 
    formula = object$formula, 
    data = newdata, 
    rules = if (object$type == "linear" || is.null(object$rules)) {NULL} else {
      structure(object$rules$description, names = object$rules$rule)}, 
    type = object$type, 
    winsfrac = winsfrac,
    x_names = object$x_names, 
    normalize = object$normalize,
    y_names = NULL,
    confirmatory = object$call$confirmatory)
  
  newdata <- newdata$x
  
  invisible(capture.output(
    rules <- print(object, penalty.par.val = penalty.par.val)))

  newdata <- newdata[, rules$rule[-which(rules$rule == "(Intercept)")]]
  if (intercept) {
    newdata <- cbind(1, newdata)
    colnames(newdata)[1] <- "(Intercept)"
  } else {
    rules <- rules[-which(rules$rule == "(Intercept)"),]
  }

  linear_terms <- rules$rule[rules$rule %in% object$x_names]
  if (length(linear_terms) > 0) {
    if (length(linear_terms) == 1L) {
      means <- mean(object$modmat[, linear_terms])
    } else {
      means <- colMeans(object$modmat[, linear_terms])  
    }
  }
  if (center.linear) {
    if (length(linear_terms) > 0) {
      newdata[, linear_terms] <- 
        t(apply(newdata[, linear_terms], 1, function(x) x - means))
    }
  }
  
  explanation <- apply(newdata, 1, function(x) x*rules$coefficient)
  rownames(explanation) <- rules$rule

  plot_func <- function(explanation, plotname, pred, maxval) {
    barplot(explanation,
            col = ifelse(explanation < 0, bar.col[1], bar.col[2]),
            main = ifelse(is.null(plotname),  
                          paste0("\n\npredicted: ", pred),
                          paste0(plotname, "\n\npredicted: ", pred)),
            xlab = ylab,
            xlim = c(-maxval, maxval), horiz = TRUE, las = 1,
            cex.axis = cex, cex.names = cex, cex.main = cex, cex.lab = cex)
  }
  maxval <- max(abs(min(explanation)), max(explanation))

  
  if (is.null(plot.obs.names)) {
    plot.obs.names <- paste("Observation", 1:nrow(newdata))
  } else if (is.logical(plot.obs.names)) {
    if (plot.obs.names) {
      plot.obs.names <- rownames(newdata)
    } else {
      plot.obs.names <- NULL
    }
  } else if (is.character(plot.obs.names)) {
    plot.obs.names <- rep_len(plot.obs.names, length.out = nrow(newdata))
  }

  par(mfrow = plot.dim)
  for (i in 1:min(nrow(newdata), plot.max.nobs)) {
    midpoints <- plot_func(explanation[ , i], plot.obs.names[i], preds[i], maxval)
    labels <- rules$description
    if (length(linear_terms) > 0) {
      linear_term_ids <- which(rules$rule %in% linear_terms)
      for (j in linear_term_ids) {
        if (center.linear) {
          labels[j] <- paste0(rules$rule[j], " = ", 
                              round(newdata[i, j], digits = digits), 
                              " (centered)")          
        } else {
          labels[j] <- paste0(rules$rule[j], " = ", 
                              round(newdata[i, j], digits = digits), 
                              " (mean = ", 
                              round(means[rules$rule[j]], digits = digits), ")")
        }
      }
    }
    if (!is.null(rule.col)) {
      text(y = midpoints, x = -maxval, labels = labels, cex = cex,
         col = rule.col, pos=4)
    }
  }
  
  return(list(predictors = newdata, contribution = explanation))
}
