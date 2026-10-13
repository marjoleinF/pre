## TODO: Fix caret method, is does not seem to use penalty.par.val properly.
## TODO: Implement (unit)  tests for pre(..., randomForest = TRUE) and get_rf_rules()
## TODO: Implement local variable importances, computed as observation level importances as per RuleFit paper.
## TODO: Implement condition filtering for rules involving multiple splits on the same variable
## TODO: Shorten code for argument checking: 
##          Group variables and check each group like: 
##          L <- list(A, B, C) all(sapply(L, class) == "matrix") or any(is.na(unlist(L)))
## TODO: Implement functionality for missing-in-attributes approach for missing values
##       MIA is supported by ctree. Can then do mean imputation before glmnet. 
##       Will need some work with predict method?
## TODO: Implement unit tests for correct passing of gamma argument
## TODO: Implement unit tests for function interact() without using nullmods
## TODO: Double check computation of H statistic


utils::globalVariables("%dopar%")

#' Derive a prediction rule ensemble
#'
#' Function \code{pre} derives a sparse ensemble of rules and/or linear functions for 
#' prediction of a continuous, binary, count, multinomial, multivariate 
#' continuous or survival response.
#' 
#' @param formula a symbolic description of the model to be fit of the form 
#' \code{y ~ x1 + x2 + ... + xn}. Response (left-hand side of the formula) 
#' should be of class numeric (for \code{family = "gaussian"} or
#' \code{"mgaussian"}), integer (for \code{family = "poisson"}), factor (for 
#' \code{family = "binomial"} or \code{"multinomial"}). See Examples below. 
#' Note that the minus sign (\code{-}) may not be used in the formula to omit
#' the intercept or variables in \code{data}, and neither should  \code{+ 0} 
#' be used to omit the intercept. To omit the intercept from the final ensemble, 
#' add \code{intercept = FALSE} to the call (although omitting the intercept from
#' the final ensemble will only very rarely be appropriate). To omit variables 
#' from the final ensemble, make sure they are excluded from \code{data}.
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
#' the program will try to make an informed guess, based on the 
#' class of the response variable specified in \code{formula}. Also see Examples 
#' below. 
#' @param ad.alpha Alpha value to be used for computing the penalty weights for the
#' adaptive lasso. Defaults to \code{NA}, yielding standard lasso estimation. 
#' To use adaptive lasso, specify a value between (and including) 0 and 1. A value of
#' 0 will yield ridge-estimated penalty weights for computing the final (lasso) 
#' penalized model. See \code{vignette("relaxed", "pre")} or 
#' argument \code{alpha} of \code{\link[glmnet]{cv.glmnet}}.
#' @param ad.penalty Penalty parameter value to be used for computing the penalty 
#' weights for the adaptive lasso. Defaults to \code{"lambda.min"}. If OLS instead
#' of elastic net regression should be used to compute weights, specify 
#' \code{ad.penalty = 0}. See also \code{vignette("relaxed", "pre")}.
#' @param use.grad logical. Should gradient boosting with regression trees be
#' employed when \code{learnrate > 0}? If \code{TRUE}, use trees fitted by 
#' \code{\link[partykit]{ctree}} or \code{\link[rpart]{rpart}} as in Friedman 
#' (2001), but without the line search. If \code{use.grad = FALSE}, 
#' \code{\link[partykit]{glmtree}} instead of \code{\link[partykit]{ctree}} 
#' will be employed for rule induction, yielding longer computation times, 
#' higher complexity, but possibly higher predictive accuracy. See Details for 
#' supported combinations of \code{family}, \code{use.grad} and \code{learnrate}.
#' @param weights optional vector of observation weights to be used for 
#' fitting the ensemble.
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
#' @param randomForest use function \code{\link[randomForest]{randomForest}}
#' of package of the same name for rule induction? If set to \code{TRUE}, 
#' arguments \code{weights}, \code{ntree}, \code{mtry} and \code{maxdepth} 
#' will be passed to
#' \code{\link[randomForest]{randomForest}}. Argument \code{learnrate} will be 
#' fixed to 0 and \code{sampfrac} will be ignored (preferences about sampling
#' procedures should be passed to \code{tree.control} when using 
#' \code{\link[randomForest]{randomForest}}) which also allows to specify other 
#' arguments of function \code{\link[randomForest]{randomForest}}.
#' @param confirmatory character vector. Specifies one or more confirmatory terms 
#' to be included in the final ensemble. Linear terms can be specified as the 
#' name of a predictor variable included in \code{data}, rules can be specified
#' as, for example, \code{"x1 > 6 & x2 <= 8"}, where x1 and x2 should be names
#' of variables in \code{data}. Terms thus specified will be included in the
#' final ensemble, as their coefficients will not be penalized in the estimation.
#' @param singleconditions \code{TRUE}, \code{FALSE} or \code{"only"}. Should 
#' rules with multiple conditions be disaggregated? Node membership for all tree
#' except the root are coded as multi-condition rules. The conditions of these
#' rules can be disaggregated to avoid selection of multi-condition rules. If
#' \code{FALSE} (the default), all non-root nodes will be included as multi-condition 
#' rules in the initial ensemble. If \code{TRUE}, all nodes will additionally be 
#' included as single-condition rules. If \code{"only"}, all nodes will be included 
#' as single-condition rules only. 
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
#' SD. \code{normalize = FALSE} will give more preference to linear terms for 
#' selection. 
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
#' @param tree.control a \code{list} with control parameters to be passed to the tree 
#' fitting function, generated using \code{\link[partykit]{ctree_control}},
#' \code{\link[partykit]{mob_control}} (if \code{use.grad = FALSE}), 
#' \code{\link[rpart]{rpart.control}} (if \code{tree.unbiased = FALSE}). Or
#' a list containing (n)one or more arguments that can be passed to function
#' \code{\link[randomForest]{randomForest}} (if \code{randomForest = TRUE}.
#' @param tree.unbiased logical. Should an unbiased tree generation algorithm 
#' be employed for rule generation? Defaults to \code{TRUE}, if set to 
#' \code{FALSE}, rules will be generated employing the CART algorithm
#' (which suffers from biased variable selection) as implemented in 
#' \code{\link[rpart]{rpart}}. See details below for possible combinations 
#' with \code{family}, \code{use.grad} and \code{learnrate}.
#' @param sparse logical. Should sparse design matrices be used? May improve
#' computation times for large datasets.
#' @param ... Further arguments to be passed to
#' \code{\link[glmnet]{cv.glmnet}}.
#' 
#' @details Note that observations with missing values will be removed prior to 
#' analysis (and a warning printed). See \code{vignette("missingness", "pre")} 
#' an \code{\link{mi_pre}} for suggestions on how to better deal with missing values.
#' 
#' In some cases, duplicated variable names may occur in the final ensemble model.
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
#' Alternatively, specifying \code{randomForest = TRUE} overrides all of these arguments and
#' employs function \code{\link[randomForest]{randomForest}} from the package of the same
#' name for rule induction. All arguments of function \code{\link[randomForest]{randomForest}}
#' should be passed using the \code{tree.control} argument, except for arguments \code{formula}, 
#' \code{data}, \code{mtry} and \code{ntrees}, because the other arguments of 
#' function \code{\link[randomForest]{randomForest}} do not have identical counterparts in 
#' function \code{\link{pre}}.
#' 
#' If an error along the lines of 'factor ... has new levels ...' is encountered, 
#' consult \code{?rare_level_sampler} for explanation and solutions.
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
#' methods for objects of class \code{pre} take a \code{penalty.parameter.val} 
#' argument, which can be used to select a different criterion.
#' 
#' If only a set of rules needs to be generated, but the final regression model
#' should not be fitted, specify the hidden argument \code{fit.final = FALSE}.
#' 
#' @examples \donttest{## Fit pre to a continuous response:
#' airq <- airquality[complete.cases(airquality), ]
#' set.seed(42)
#' airq.ens <- pre(Ozone ~ ., data = airq)
#' airq.ens
#' ## Use relaxed lasso to estimate final model
#' airq.ens.rel <- pre(Ozone ~ ., data = airq, relax = TRUE)
#' airq.ens.rel
#' ## Use adaptive lasso to estimate final model
#' airq.ens.ad <- pre(Ozone ~ ., data = airq, ad.alpha = 0)
#' airq.ens.ad
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
#' iris.ens <- pre(Species ~ ., data = iris, family = "multinomial")
#' iris.ens
#' 
#' ## Fit pre to a survival response:
#' library("survival")
#' lung <- lung[complete.cases(lung), ]
#' set.seed(42)
#' lung.ens <- pre(Surv(time, status) ~ ., data = lung, family = "cox")
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
#' @seealso \code{\link{pre}}, \code{\link{importance.pre}}
#' @references Fokkema, M. (2020). Fitting prediction rule ensembles with R 
#' package pre. \emph{Journal of Statistical Software, 92}(12), 1-30.
#' \doi{10.18637/jss.v092.i12}
#' 
#' Fokkema, M. & Strobl, C. (2020). Fitting prediction rule ensembles to psychological 
#' research data: An introduction and tutorial. \emph{Psychological Methods 25}(5), 
#' 636-652. \doi{10.1037/met0000256}, \url{https://arxiv.org/abs/1907.05302}
#' 
#' Friedman, J. H. (2001). Greedy function approximation: a gradient boosting 
#' machine. \emph{The Annals of Applied Statistics, 29}(5), 1189-1232.
#' 
#' Friedman, J. H., & Popescu, B. E. (2008). Predictive learning via rule 
#' ensembles. \emph{The Annals of Applied Statistics, 2}(3), 916-954, \doi{10.1214/07-AOAS148}.
#' 
#' Hothorn, T., & Zeileis, A. (2015). partykit: A modular toolkit for recursive 
#' partytioning in R. \emph{Journal of Machine Learning Research, 16}, 3905-3909.
#' 
#' @export
pre <- function(formula, data, family = gaussian, ad.alpha = NA, ad.penalty = "lambda.min",
                use.grad = TRUE, weights, type = "both", sampfrac = .5, 
                maxdepth = 3L, learnrate = .01, mtry = Inf, ntrees = 500,
                confirmatory = NULL, singleconditions = FALSE,
                randomForest = FALSE, 
                winsfrac = .025, normalize = TRUE, standardize = FALSE,
                ordinal = TRUE, nfolds = 10L, tree.control, tree.unbiased = TRUE, 
                removecomplements = TRUE, removeduplicates = TRUE, 
                verbose = FALSE, par.init = FALSE, par.final = FALSE, 
                sparse = FALSE, ...) { 
  
  
  #####################
  ## Check arguments ##
  #####################
  
  ## Save call:
  cl <- match.call()
  
  ## Check if proper formula argument is specified 
  if (!(inherits(formula, "formula"))) {
    stop("Argument formula should specify and object of class 'formula'.\n")
  }
  ## Check if dot and functions are simultaneously used in formula
  form <- as.character(formula[3])
  for (i in names(data)) {
    form <- gsub(pattern = i, replacement = "", x = form)
  }
  if (any(grepl(".", form, fixed = TRUE))) {
    if (any(grepl("(", form, fixed = TRUE))) {
      if (any(grepl(")", form, fixed = TRUE))) {
        warning("Argument formula contains both one or more functions of predictor variables, as well as a dot ('.'), which should be avoided. Model fitting may fail, and/or both the original variable(s) and their functions may be included as predictor variables.\n", immediate. = TRUE)  
      }
    }
  }
  if (any(grepl("-", as.character(formula), fixed = TRUE))) {
    warning("Argument formula contains a minus sign. Note that the minus sign should not be used to omit the intercept or variables from the ensemble. To omit the intercept from the final ensemble, specify intercept = FALSE\n", immediate. = TRUE)
  }

  ## Check if proper data argument is specified:
  if (!is.data.frame(data)) {
    stop("Argument data should specify a data frame.\n")
  }

  ## Check and set up family argument: 
  if (length(family) > 1L) {
    warning("Argument family has length > 1, only first element will be used.\n")
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
      stop("Argument family should specify a family object.\n")
    }
  }
  if (is.character(family)) {
    if (!(family %in% c("gaussian", "binomial", "poisson", "mgaussian", "multinomial", "cox"))) {
      stop("Argument family should be equal to 'gaussian', 'binomial', 'poisson', 'multinomial', 'mgaussian', 'cox', or a corresponding family object.\n")
    }
  } else {
    stop("Argument family should be equal to 'gaussian', 'binomial', 'poisson', 'multinomial', 'mgaussian', 'cox', or a corresponding family object.\n")
  }
  
  ## Check if proper weights argument is specified, if specified:
  if (missing(weights)) {
    weights <- rep(1L, times = nrow(data))
  } else {
    if (!is.numeric(weights)) {
      warning("Argument weights should, but does not specify a vector of numeric values.\n", immediate. = TRUE)
    }
    if (length(weights) != nrow(data)) {
      warning("Length of argument weights should be, but is not equal to nrow(data).\n", immediate. = TRUE)
    }
    if (!all(weights > 0)) {
      warning("All weights should be positive but some are not.\n", immediate. = TRUE)
    }
  }
  
  ## Check if proper type argument is specified:
  if (!(length(type) == 1L && type %in% c("rules", "both", "linear"))) {
    stop("Argument type should be 'rules', 'linear' or 'both'.\n")
  }
  
  ## Check if user is trying to specify a random forest as in example Fokkema 2020 JSS
  if (mtry < Inf && tree.unbiased == FALSE) {
    message("Use of mtry argument with tree.unbiased = FALSE is not supported. Alternatively, it is supported through Breiman's (2000) original random forest algorithm, which can be used for rule generation by specifying randomForest = TRUE. See ?pre and argument randomForest for further info.")
  }

  
  ## Check if proper sampfrac argument is specified:
  if (!is.function(sampfrac)) {
    if (!(length(sampfrac) == 1 && is.numeric(sampfrac) && sampfrac > 0 && 
          sampfrac <= 1)) {
      stop("Argument sampfrac should be a single numeric value > 0 and <= 1, or a sampling function.\n")
    }
  }

  ## Check if proper maxdepth argument is specified:
  if (is.function(maxdepth)) {
    maxdepth <- maxdepth(ntrees = ntrees)
  } else if (!is.numeric(maxdepth)) {
    stop("Argument maxdepth should be either a numeric vector of length 1 or ntrees, or a random number generating function.\n")
  } else if (!(length(maxdepth) %in% c(1L, ntrees))) {
    warning("Argument maxdepth should be either a numeric vector of length 1 or ntrees, only first value of maxdepth will be used.\n")
    maxdepth <- maxdepth[1L]
  } 
  if (!all(maxdepth > 0)) {
    stop("All values of maxdepth should be > 0.\n")
  } 
  if (!all(maxdepth == suppressWarnings(as.integer(maxdepth)) | is.infinite(maxdepth))) {
    stop("Argument maxdepth should consist of  of integer values or Inf), or a random number generating function.\n")
  }
  
  ## Check if proper learnrate argument is specified:
  if (!(length(learnrate) == 1L && is.numeric(learnrate) && 
        (learnrate >= 0 || learnrate <= 1))) {
    stop("Argument learnrate shoud be a single numeric value >= 0 and <= 1.\n")
  }
  
  ## Check if proper mtry argument is specified:
  if (!(length(mtry) == 1L && mtry > 0 && 
        (mtry == suppressWarnings(as.integer(mtry)) || is.infinite(mtry)))) {
    stop("Argument mtry should be a single integer value, or Inf.\n")
  }
  
  ## Check if proper ntrees argument is specified:
  if (!(length(ntrees) == 1L && ntrees == as.integer(ntrees) && ntrees > 0)) {
    stop("Argument ntrees should be a single positive integer.\n")
  }
  
  ## Check if proper winsfrac argument is specified:
  if (!(length(winsfrac == 1L) && is.numeric(winsfrac) && winsfrac >= 0 && 
        winsfrac < 1)) {
    stop("Argument winsfrac should be a numeric value >= 0 and < 1.\n")
  }

  ## Check if proper nfolds argument is specified:
  if (!(length(nfolds) == 1L && nfolds > 0 && nfolds == as.integer(nfolds))) {
    stop("Argument nfolds should be a positive integer.\n")
  }
  
  ## Check if proper confirmatory argument is specified:
  if (!is.null(confirmatory)) {
    if (!is.character(confirmatory)) {
      stop("Argument confirmatory should specify a character vector.\n")
    }
  }
  
  ## Check if logical arguments of length 1 are properly specified:
  is_logical_and_length_one <- function(x) {is.logical(x) && length(x) == 1L}
  for (i in c(use.grad, removeduplicates, removecomplements, normalize, 
              standardize, ordinal, verbose, tree.unbiased, par.init, 
              par.final)) {
    if (!is_logical_and_length_one(i)) {
      stop("Argument ", i, "should be TRUE or FALSE.\n")
    }
  }
  
  if (par.final || par.init) {
    if (!requireNamespace("foreach")) {
      warning("Parallel computation requires package foreach. Arguments par.init and par.final are set to FALSE.\n")   
      par.init <- par.final <- FALSE
    }
  }

  ## Check if proper tree.control argument is specified:
  if (missing(tree.control)) {
    if (randomForest) {
      tree.control <- list(mtry = mtry, ntree = ntrees)
    } else {
      if (tree.unbiased && use.grad) {
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
          warning("Value specified for mtry will be ignored if tree.unbiased = FALSE and randomForest = FALSE.\n")
        }
      }
    } 
  } else {
    if (!is.list(tree.control)) {
      stop("Argument tree.control should specify a list of control parameters.\n")
    }
    if (randomForest) {
      if (is.null(tree.control$mtry)) tree.control$mtry <- mtry
      if (is.null(tree.control$ntree)) tree.control$ntree <- ntrees
    } else if (use.grad && tree.unbiased) {
      if (!setequal(names(ctree_control()), names(tree.control))) {
        args_not_in_ctree_control <- names(tree.control)[which(
          !names(tree.control) %in% names(ctree_control()))]
        warning("Argument tree.control contains named elements: ", 
                paste0(args_not_in_ctree_control, collapse = ", "), 
                ", which should probably not be there, see ?ctree_control.\n")
      }
    } else if (!use.grad && tree.unbiased) { 
      if (!setequal(names(mob_control()), names(tree.control))) {
        args_not_in_glmtree_control <- names(tree.control)[which(
          !names(tree.control) %in% names(mob_control()))]
        warning("Argument tree.control has named elements ", 
                 paste0(args_not_in_glmtree_control, collapse = ", "), 
                ", which should probably not be there, see ?mob_control.\n")
      }
    } else if (!tree.unbiased) {
      if (!setequal(names(rpart.control()), names(tree.control))) {
        warning("Argument tree.control should be a list containing names elements ", paste(names(rpart.control()), collapse = ', '), "\n")
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
      warning("Argument tree.unbiased was set to FALSE, so rpart is employed for tree induction, and value specified for mtry will be ignored.\n")
      mtry <- Inf
    }
  }
  
  if (!tree.unbiased) {
    if (!requireNamespace("rpart")) {
      warning("Arguments tree.unbiased is set to FALSE, but package rpart is not installed. An error will likely occur next. Install rpart first, then rerun.\n")   
    }
    if (!use.grad && learnrate > 0) {
      stop("Employing the rpart algorithm with a learnrate > 0 without gradient boosting is not supported.\n")
    }
  }
  
  
  ######################################
  ## Prepare data, formula and family ##
  ######################################

  if (!is.null(tree.control$cluster)) {
    ## If cluster variable has been specified, make sure it is included in data 
    if (!inherits(tree.control$cluster, "name")) 
      warning("A cluster argument was passed to argument tree.control, but it is not a name, so it might be ignored.")
    cluster <- data[ , as.character(tree.control$cluster)]
  }
  
  ## prepare model frame:
  data <- model.frame(Formula::as.Formula(formula), data = data, na.action = NULL)
  if (!is.null(tree.control$cluster)) {
    ## If cluster variable has been specified, make sure it is included in data 
    data[ , as.character(tree.control$cluster)] <- cluster
  }

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
  
  ## get response variable name(s):
  y_names <- names(data)[attr(attr(data, "terms"), "response")]
  if (family == "mgaussian" || length(y_names) == 0) {
    y_names <- attr(terms(Formula(formula), rhs = 0, data = data), "term.labels")
    family <- "mgaussian"
  }
  
  ## get predictor variable names:
  if (family == "cox" || is.Surv(data[ , y_names])) {
    x_names <- attr(attr(data, "terms"), "term.labels")
  } else {
    x_names <- attr(terms(Formula(formula), lhs = 0, data = data), "term.labels")
  }
  
  ## Coerce ordered categorical variables to numeric:
  if (ordinal) {
    if (any(ord_var_inds <- sapply(data[ , x_names], is.ordered))) {
      data[ , names(which(ord_var_inds))] <- sapply(
        data[ , names(which(ord_var_inds))], as.numeric)
    }
  }
  
  ## expand dot and put ticks around variables within functions, if present:
  if (family != "mgaussian") {
    formula <- formula(data)
  } else {
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
      if (is.factor(data[ , y_names])) { # then family should be bi- or multinomial
        if (is.ordered(data[ , y_names])) {
          warning("An ordered factor was specified as the response variable, which will be treated as an unordered factor response.")
          data[ , y_names] <- factor(data[ , y_names], ordered = FALSE)
        } 
        if (nlevels(data[ , y_names]) == 2L) {
            family <- "binomial"
        } else if (nlevels(data[ , y_names]) > 2L) {
            family <- "multinomial"
          }
      } else if (is.Surv(data[ , y_names])) { # then family should be cox
        family <- "cox"
      } else if (!is.numeric(data[ , y_names])) { # then response is not a factor, survival or numeric
        warning("The response variable specified through argument formula should be of class numeric, factor or Surv.")
      }
    } else if (length(y_names) > 1L) { # multiple responses specified, should be numeric
      if (all(sapply(data[ , y_names], is.numeric))) {
        family <- "mgaussian"
      } else {
        warning("Multiple response variables were specified, but not all were (but should be) numeric.\n")
      }
    }
    
  } else { # family was specified, check if correct;
    
    if (family[1L] == "gaussian") {
      if (length(y_names) > 1L) {
        warning("Argument family was set to 'gaussian', but multiple response variables were specified in formula. Consider specifying family = 'mggaussian'?\n")        
      }
      if (!is.numeric(data[ , y_names])) { # then family should be poisson or gaussian
        warning("Argument family was set to 'gaussian', but the response variable specified in formula is not of class numeric.\n")
      }
    } else if (family[1L] == "poisson") {
      if (length(y_names) > 1L) {
        warning("Argument family was set to 'poisson', but multiple response variables were specified, which is not supported.\n")
      }
      if (!isTRUE(all.equal(round(data[ , y_names]), data[ , y_names]))) {
        warning("Argument family' was set to 'poisson', but the response variable specified in formula is non-integer.\n")
      }
    } else if (family[1L] == "multinomial") {
      if (length(y_names) > 1L) {
        warning("Argument family was set to 'multinomial', but multiple response variables were specified in formula, which is not supported. Check specified response variable (should be a single factor with > 2 levels) and family.\n")        
      }
      if (!is.factor(data[ , y_names])) {
        warning("Argument family was set to 'multinomial', but response variable is numeric. Response variable will be converted to factor.")
        data[ , y_names] <- factor(data[ , y_names])
      }  
    } else if (family[1L] == "binomial") {
      if (length(y_names) > 1L) {
        warning("Argument family was set to 'binomial', but multiple response variables were specified, which is not supported.\n")
      } else if (!is.factor(data[ , y_names])) {
        warning("Argument family was set to 'binomial', but the response variable specified is not a factor.\n")
      } else if (is.ordered(data[ , y_names])) {
        warning("Argument family was set to 'binomial', but the response variable specified is an ordered factor. It will be treated as an unordered factor.")
      } else if (nlevels(data[ , y_names]) != 2L) {
        warning("Argument family was set to 'binomial', but the response variable has ", nlevels(data[ , y_names]), " levels.\n")        
      }
    } else if (family[1L] == "multinomial") {
      if (length(y_names) > 1L) {
        warning("Argument family was set to 'multinomial', but multiple response variables were specified, which is not supported.\n")
      } else if (!is.factor(data[ , y_names])) {
        warning("Argument family was set to 'multinomial', but the response variable specified is not a factor.\n")
      } else if (is.ordered(data[ , y_names])) {
        warning("Argument family was set to 'multinomial', but the response variable specified is an ordered factor. It will be treated as an unordered factor.")
      } else if (nlevels(data[ , y_names]) < 3L) {
        warning("Argument family was set to 'multinomial', but the response variable has ", nlevels(data[ , y_names]), " levels.\n")
      }
    } else if (family[1L] == "cox") {
      if (length(y_names) > 1L) {
        warning("Argument family was set to 'cox', but multiple response variables were specified, which is not supported.\n")
      } else if (!is.Surv(data[ , y_names])) {
        warning("Argument family was set to 'cox', but the response variable specified is not of class Surv.\n")
      }
    } else if (family == "mgaussian") {
      if (length(y_names) == 1L) {
        warning("Argument family was set to 'mgaussian', but only a single response variable was specified.\n")
      } else if (!all(sapply(data[ , y_names], is.numeric))) {
        warning("Argument family was set to 'mgaussian', but not all response variables specified are numeric.\n")
      }
    }
  }

  ## Check specification of tree growing algorithms employed:
  if (!tree.unbiased) { # rpart is employed
    if (family == "mgaussian") {
      stop("Employing rpart algorithm for rule induction with a multivariate response variable is not supported. Specify tree.unbiased = TRUE and use.grad = FALSE.\n")
    } else if (learnrate > 0 && family == "multinomial") {
      stop("Employing rpart algorithm for rule induction with a multinomial response variable and learnrate > 0 is not supported. Specify learnrate = 1, or tree.unbiased = TRUE and use.grad = TRUE.\n")
    }
  } else if (!use.grad) { # (g)lmtree is employed
    if (family == "multinomial") {
      stop("Employing (g)lmtree for rule induction with a multinomial response variable is not supported. Specify use.grad = TRUE for multivariate responses.\n")
    } else if (family == "mgaussian") {
      stop("Employing (g)lmtree for rule induction with a multivariate response variable is not supported. Specify use.grad = TRUE for multivariate responses.\n")
    } else if (family == "cox") {
      stop("Employing (g)lmertree for rule induction with a survival response is not supported. Specify use.grad = TRUE for a survival response.\n")
    }
  }
  
  if (randomForest && !family %in% c("gaussian", "multinomial", "binomial")) 
    stop("randomForest can only be used for rule generation for numeric and categorical response variables (i.e., gaussian, binomial or multinomial family).") 
    
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
    if (isTRUE(all.equal(round(data[ ,  y_names]), data[ ,  y_names]))) { # if response passes integer test
      data[ ,  y_names] <- data[ ,  y_names] + 0.01 # add small constant to response to prevent response being interpreted as count by ctree or rpart
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
      cat("\nA rule ensemble for prediction of a continuous response will be created.\n")
    } else if (family == "poisson") {     
      cat("\nA rule ensemble for prediction of a count response will be created.\n")
    } else if (family == "binomial") {
      cat("\nA rule ensemble for prediction of a binary categorical response will be created.\n")
    } else if (family == "multinomial") {
      cat("\nA rule ensemble for prediction of a multinomial response will be created.\n")
    } else if (family == "mgaussian") {
      cat("\nA rule ensemble for prediction of a multivariate continuous response will be created.\n")
    } else if (family == "cox") {
      cat("\nA rule ensemble for prediction of a survival response will be created.\n")
    }
  }
  
    
  #############################
  ## Derive prediction rules ##
  #############################
  
  if (type == "linear") {
    rules <- NULL
    rulevars <- NULL
  } else {
    list_dots <- list(...)
    offset <- if (is.null(list_dots)) NULL else list_dots$offset
    rule_object <- try(pre_rules(formula = formula, 
                                 data = data,
                                 offset = offset,
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
                                 singleconditions = singleconditions,
                                 removeduplicates = removeduplicates, 
                                 removecomplements = removecomplements,
                                 tree.unbiased = tree.unbiased,
                                 return.dupl.compl = TRUE, 
                                 sparse = sparse,
                                 randomForest = randomForest))
    if (inherits(rule_object, "try-error")) {
      if (grepl("has new levels", rule_object)) {
        stop(rule_object[1], paste("\n Hint: There may be a predictor variable which is a factor with rare levels, this likely yields problems for rule estimation. Consult ?rare_level_sampler"))
      }
      stop(rule_object[1])
    }
    rules <- rule_object$rules
    rulevars <- rule_object$rulevars
  }

  #################################################################
  ## Prepare rules, linear terms, outcome variable for selection ##
  #################################################################
  
  if (type == "rules" && length(rules) == 0) {
    warning("No prediction rules could be derived from the data, yielding an intercept-only model. Estimation will be halted.")
    return(NULL)
  }
  
  if (is.numeric(small_constant_added)) {
    data[ , y_names] <- data[ , y_names] - small_constant_added
  }
  
  modmat_data <- get_modmat(
    formula = formula, 
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
  
  y <- modmat_data$y
  x <- modmat_data$x  

  if (is.null(cl$penalty.factor)) {
    penalty.factor <- rep(1L, times = ncol(x))
  } 
  if (!is.null(confirmatory)) {
    penalty.factor[which(colnames(x) %in% confirmatory)] <- 0L
  }

  ## check for duplicates in variable names:
  ## (can happen, for example, due to labeling of dummy indicators for factors)
  if (!(length(unique(colnames(x))) == length(colnames(x)))) { 
    warning("There are variables in the model with overlapping variable names. If predictor variables of type factor were specified with numbers in their name, consider renaming these and rerunning the analysis. See 'Details' under ?pre.") 
  } 
  
  if (!is.null(cl$fit.final)) { ## Hidden argument has been used
    if (!cl$fit.final) {
      glmnet.fit <- NULL 
      if (verbose) cat("\n\nThe final penalized regression model will not be fit.\n")
    } else {
      warning("Hidden argument fit.final was specified, but it was not set to FALSE and will likely be ignored.")
    }
  } else {
    if (is.na(ad.alpha)) {
      glmnet.fit <- cv.glmnet(x, y, nfolds = nfolds, weights = weights, 
                              family = family, parallel = par.final, 
                              standardize = standardize, 
                              penalty.factor = penalty.factor,
                              ...)
    } else { ## use adaptive lasso, apply cv.glmnet twice
      
      cv.glmnet.args <- list(...)
      ## Fix fold ids for reuse
      cv.glmnet.args$foldid <- if (is.null(cl$foldid)) 
        sample(rep(1:nfolds, length = nrow(data))) else cl$foldid
      cv.glmnet.args$x <- x
      cv.glmnet.args$y <- y
      cv.glmnet.args$weights <- weights
      cv.glmnet.args$family <- family
      cv.glmnet.args$parallel <- par.final
      cv.glmnet.args$standardize <- standardize
      cv.glmnet.args$alpha <- ad.alpha
      cv.glmnet.args$penalty.factor <- penalty.factor
      cv.glmnet.args$relax <- FALSE
      
      ## Compute penalty weights
      glmnet.fit <- do.call(cv.glmnet, cv.glmnet.args)
      penalty.factor <- coef(glmnet.fit, s = ad.penalty)
      if (is.list(penalty.factor)) {
        penalty.factor <- rowMeans(do.call(cbind, penalty.factor))
        penalty.factor <- 1L / abs(as.numeric(penalty.factor)[
          -which(names(penalty.factor) == "(Intercept)")])
      } else {
        penalty.factor <- 1L / abs(as.numeric(penalty.factor)[
          -which(rownames(penalty.factor) == "(Intercept)")])
      }
      if (!is.null(cl$confirmatory)) {
        penalty.factor[which(colnames(x) %in% cl$confirmatory)] <- 0L
      }
      
      ## Fit final model
      cv.glmnet.args$alpha <- if (is.null(cl$alpha)) 1 else cl$alpha
      cv.glmnet.args$penalty.factor <- penalty.factor
      cv.glmnet.args$relax <- if (is.null(cl$relax)) FALSE else cl$relax
      glmnet.fit <- do.call(cv.glmnet, cv.glmnet.args)

    }
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
  }
  
  
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