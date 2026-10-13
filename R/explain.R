



#' Explain predictions from final prediction rule ensemble
#'
#' \code{explain} shows which rules apply to which observations and visualizes
#' the contribution of rules and linear predictors to the predicted values
#'
#' @param object object of class \code{\link{pre}}.
#' @param newdata optional dataframe of new (test) observations, including all
#' predictor variables used for deriving the prediction rule ensemble.
#' @inheritParams print.pre
#' @param response numeric or character vector of length one. Specifies the
#' name or number of the response variable (for multivariate responses) or
#' the name or number of the factor level (for multinomial responses) for 
#' which explanations and contributions should be computed and/or plotted.
#' Only used for\code{pre}s fitted to multivariate or multinomial responses.  
#' @param plot logical. Should explanations be plotted?
#' @param intercept logical. Specifies whether intercept should be included in
#' explaining predictions.
#' @param center.linear logical. Specifies whether linear terms should be
#' centered with respect to the training sample mean before computing their 
#' contribution to the predicted value. If \code{intercept = TRUE}, this
#' will also affect the intercept. That is, the value of the intercept returned
#' will differ from that of the value returned by the \code{print} method.
#' @param pred.type character. Specifies the type of predicted values to be 
#' computed, returned and provided in the plot(s). Note that the computed 
#' contributions must be additive and are therefore always on the scale of 
#' the linear predictor.
#' @param plot.max.nobs numeric. Specifies maximum number of observations
#' for which explanations will be plotted. The default (\code{4}) plots the
#' explanation for the first four observations supplied in \code{newdata}.
#' @param plot.dim numeric vector of length 2. Specifies the number of rows and
#' columns in the resulting plot.
#' @param cex numeric. Specifies the relative text size of title, tick and axis
#' labels.
#' @param deprecation.message logical. Should deprecation message be returned?
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
#' @param ... Further arguments to be passed to \code{\link{predict.pre}} and 
#' \code{\link[glmnet]{predict.cv.glmnet}}.
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
#' @references Fokkema, M. & Strobl, C. (2020). Fitting prediction rule 
#' ensembles to psychological research data: An introduction and tutorial. 
#' \emph{Psychological Methods 25}(5), 636-652. \doi{10.1037/met0000256},
#' \url{https://arxiv.org/abs/1907.05302}
#' @examples \donttest{airq <- airquality[complete.cases(airquality), ]
#' set.seed(1)
#' train <- sample(1:nrow(airq), size = 100)
#' set.seed(42)
#' airq.ens <- pre(Ozone ~ ., data = airq[train,])
#' airq.ens.exp <- explain(airq.ens, newdata = airq[-train,])
#' airq.ens.exp$predictors
#' airq.ens.exp$contribution
#' 
#' ## Can also include intercept in explanation:
#' airq.ens.exp <- explain(airq.ens, newdata = airq[-train,])
#' 
#' ## Fit PRE with linear terms only to illustrate effect of center.linear:
#' set.seed(42)
#' airq.ens2 <- pre(Ozone ~ ., data = airq[train,], type = "linear")
#' ## When not centered around their means, Month has negative and 
#' ##   Day has positive contribution:
#' explain(airq.ens2, newdata = airq[-train,][1:2,],
#'         penalty.par.val = "lambda.min")$contribution
#' ## After mean centering, contributions of Month and Day have switched
#' ##   sign (for these two observations): 
#' explain(airq.ens2, newdata = airq[-train,][1:2,], 
#'         penalty.par.val = "lambda.min", center.linear = TRUE)$contribution
#' }
#' @seealso \code{\link{shap}}, \code{\link{pre}}, \code{\link{plot.pre}},
#' \code{\link{coef.pre}}, \code{\link{importance.pre}}, \code{\link{cvpre}},
#' \code{\link{interact}}, \code{\link{print.pre}}
#' @export
explain <- function(object, newdata, penalty.par.val = "lambda.1se",
                    response = 1L, plot = TRUE, intercept = FALSE, 
                    center.linear = FALSE, plot.max.nobs = 4, 
                    plot.dim = c(2, 2), plot.obs.names = TRUE, 
                    pred.type = "response", digits = 3L,  cex = .8,
                    ylab = "Contribution to linear predictor",
                    bar.col = c("#E495A5", "#39BEB1"),
                    rule.col = "darkgrey", deprecation.message = TRUE, ...) {
  
  if (deprecation.message) {
    .Deprecated(
      new = "shap",
      msg = "Function explain() is deprecated. Function shap() is now available since version 1.10.0 for explaining predictions through SHAP values. Users are adviced to compute SHAP values as these have much more favorable properties than the results returned by function explain()."
  )}
  
  ## check arguments
  if (!inherits(object, what = "pre")) {
    stop("Argument object should specify an object of class 'pre'.")
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
  
  ## Check if all variables with non-zero importance have been supplied:
  req_vars <- importance(object, plot = FALSE, 
                         penalty.par.val = penalty.par.val)$varimps$varname
  if (!all(req_vars%in% names(newdata))) {
    stop("All variables with non-zero importances should be included in newdata")
  }
  
  ## Check if proper response variable is specified
  if (object$family %in% c("mgaussian", "multinomial")) {
    if (!((is.numeric(response) || is.character(response)) && length(response) == 1L)) {
      stop("Argument response should specify a numeric or character vector or length 1.")
    }
    if (object$family == "mgaussian") {
      if (is.numeric(response)) {
        if (response > length(object$y_names)) {
          stop(paste0("There is no response variable number ", response, "."))
        }
        response <- object$y_names[response]
      } else {
        if (!(response %in% object$y_names)) {
          stop(paste0("Response variable named '", response, "'does not exist. Argument response should specify one of: "), paste(object$y_names, collapse = " "))
        }
      }
    } else {
      if (is.numeric(response)) {
        if (response > length(levels(object$data[ , object$y_names]))) {
          stop(paste0("The response variable has no level number ", response, "."))
        }
        response <- levels(object$data[ , object$y_names])[response]
      } else {
        if (!(response %in% levels(object$data[ , object$y_names]))) {
          stop(paste0("Response variable level '", response, "' does not exist. Argument response should specify one of: "), paste(levels(object$data[ , object$y_names]), collapse = " "))
        }
      }
    }
  }
  
  ## Prepare model matrix for getting explanations
  modmat <- newdata
  ## Add values of variables which have zero importance, if missing:
  nonreq_vars <- object$x_names[!object$x_names %in% req_vars]
  if (any(add_vars <- !nonreq_vars %in% names(modmat))) {
    modmat[, nonreq_vars[add_vars]] <- object$data[1, nonreq_vars[add_vars]]
  }
  
  ## Get predicted values
  preds <- round(predict(object, newdata = modmat, type = pred.type,
                         penalty.par.val = penalty.par.val, ...), 
                 digits = digits)
  
  ## Prepare newdata
  winsfrac <- (object$call)$winsfrac
  if (is.null(winsfrac)) winsfrac <- formals(pre)$winsfrac
  
  ## Check if variable names and classes are the same in newdata as in object$data
  if (!all(object$x_names %in% names(modmat))) {
    modmat <- model.frame(as.Formula((object$call)$formula), data = modmat, 
                          rhs = NULL, lhs = 0, na.action = NULL)
  } else {
    modmat <- modmat[ , object$x_names]
  }
  
  ## Coerce character and logical variables to factors
  if (any(char_names <- sapply(modmat, is.character))) {
    char_names <- names(modmat)[char_names]
    data[ , char_names] <- sapply(modmat[ , char_names], factor)
  }
  if (any(logic_names <- sapply(modmat, is.logical))) {
    logic_names <- names(modmat)[logic_names]
    modmat[ , logic_names] <- sapply(modmat[ , logic_names], factor)
  } 
  
  ## Coerce ordered categorical variables to numeric, if necessary
  if (if (is.null((object$call)$ordinal)) {
    formals(pre)$ordinal
  } else {
    (object$call)$ordinal
  }) {
    if (any(ordered_names <- sapply(modmat, is.ordered))) {
      ordered_names <- names(modmat)[ordered_names]
      modmat[ , ordered_names] <- sapply(modmat[ , ordered_names], as.numeric)
    }
  }
  
  if (any(is.na(modmat))) {
    modmat <- modmat[complete.cases(modmat),]
    warning("Some observations in newdata have missing predictor variable values and will be removed.", immediate. = TRUE)
  }
  
  ## Check and set factor levels of newdata to variable levels in object$data
  if (any(factor_inds <- sapply(modmat, is.factor))) {
    for (i in names(modmat)[factor_inds]) {
      if (all(levels(modmat[ , i]) %in% levels(object$data[ , i]))) {
        levels(modmat[ , i]) <- levels(object$data[ , i])
      } else {
        stop("Variable ", i, " has levels not present in training data. Cannot compute predictions.")
      }
    }
  }
  
  modmat <- get_modmat(
    wins_points = object$wins_points, 
    x_scales = object$x_scales, 
    formula = object$formula, 
    data = modmat, 
    rules = if (object$type == "linear" || is.null(object$rules)) {NULL} else {
      structure(object$rules$description, names = object$rules$rule)}, 
    type = object$type, 
    winsfrac = winsfrac,
    x_names = object$x_names, 
    normalize = object$normalize,
    y_names = NULL,
    confirmatory = object$call$confirmatory)$x
  if (intercept) modmat <- cbind(`(Intercept)` = 1L, modmat)
  coefs <- coef(object$glmnet.fit, s = penalty.par.val)
  
  ## Select only non-zero terms and compute explanations:
  if (object$family %in% c("multinomial", "mgaussian")) {
    coefs <- sapply(coefs, function(x) x[x@i + 1L, ])[ , response]
  } else {
    coefs <- coefs[coefs@i + 1L, ]
  }
  
  if (!intercept) coefs <- coefs[-1L]
  modmat <- modmat[ , names(coefs)]      
  linear_terms <- names(coefs)[!grepl("rule", names(coefs))]
  linear_terms <- linear_terms[!grepl("(Intercept)", linear_terms)]
  numeric_linear_terms <- linear_terms[linear_terms %in% object$x_names]
  
  if (length(numeric_linear_terms) > 0L) {
    tmp <- scale(modmat[ , numeric_linear_terms], center = TRUE, scale = FALSE)
    means <- attr(tmp, "scaled:center")
    if (center.linear) modmat[ , numeric_linear_terms] <- tmp
    if (is.null(names(means))) names(means) <- numeric_linear_terms
    if (object$normalize) means <- means * object$x_scales[numeric_linear_terms]
  }
  explanation <- apply(modmat, 1L, function(x) x*coefs)
  
  ## Combine dummy indicators of categorical predictors
  factor_names <- c()
  for (i in linear_terms[!(linear_terms %in% numeric_linear_terms) & 
                         !(linear_terms %in% object$call$confirmatory)]) {
    factor_names <- c(factor_names, 
                      unlist(sapply(object$x_names, function(x) grep(x, i))))
  }
  factor_names <- unique(names(factor_names))
  ## TODO: Make this optional through an argument of function explain()?
  ## Replace all dummy indicators from the same factor_names with a single sum
  for (i in factor_names) {
    factor_ids <- grep(i, rownames(explanation))
    if (length(factor_ids) > 1L) {
      explanation[factor_ids[1L], ] <- colSums(
        explanation[grep(i, rownames(explanation)),])
      explanation <- explanation[-(factor_ids[-1L]), ]
    } else {
      explanation[factor_ids[1L], ] <- sum(
        explanation[grep(i, rownames(explanation)),])
    }
    rownames(explanation)[factor_ids[1L]] <- i
  }
  
  if (plot) {
    if (length(numeric_linear_terms) > 0L && intercept && center.linear) {
      explanation["(Intercept)", ] <- explanation["(Intercept)", ] + 
        sum((means[numeric_linear_terms] / object$x_scales[numeric_linear_terms]) * 
              coefs[numeric_linear_terms])
    }
    plot_func <- function(explanation, plotname, pred, maxval) {
      barplot(explanation,
              col = ifelse(explanation < 0L, bar.col[1L], bar.col[2L]),
              main = ifelse(is.null(plotname),  
                            paste0("\n\npredicted value: ", pred),
                            paste0(plotname, "\n\npredicted value: ", pred)),
              xlab = ylab, xlim = maxval, horiz = TRUE, las = 1L,
              cex.axis = cex, cex.names = cex, cex.main = cex, cex.lab = cex)
    }
    maxval <- c(min(explanation), max(explanation))
    
    if (is.null(plot.obs.names)) {
      plot.obs.names <- paste("Observation", 1L:nrow(newdata))
    } else if (is.logical(plot.obs.names)) {
      if (plot.obs.names) {
        plot.obs.names <- rownames(newdata)
      } else {
        plot.obs.names <- NULL
      }
    } else if (is.character(plot.obs.names)) {
      plot.obs.names <- rep_len(plot.obs.names, length.out = nrow(newdata))
    }
    
    labels <- coef(object, penalty.par.val = penalty.par.val)
    rownames(labels) <- labels$rule
    labels <- labels[rownames(explanation) , "description"]
    if (length(factor_names) > 0L) {
      if (all(rownames(explanation)[is.na(labels)] %in% factor_names)) {
        labels[is.na(labels)] <- rownames(explanation)[is.na(labels)]
      } else {
        warning("Something may have gone wrong with the computing and coding of the following variables: ", 
                paste(rownames(explanation)[!(rownames(explanation)[is.na(labels)] %in% factor_names)]))
      }
    }
    
    for (i in 1L:min(nrow(newdata), plot.max.nobs)) {
      midpoints <- plot_func(explanation[ , i], plot.obs.names[i], preds[i], maxval)
      if (length(numeric_linear_terms) > 0) {
        for (j in numeric_linear_terms) {
          label_id <- which(rownames(explanation) == j)
          if (center.linear) {
            labels[label_id] <- paste0(j, " = ", round(newdata[i, j], digits = digits), 
                                       " (centered)")          
          } else {
            labels[label_id] <- paste0(j, " = ", round(newdata[i, j], digits = digits), 
                                       " (mean = ", 
                                       round(means[j], digits = digits), ")")
          }
        }
        for (j in factor_names) {
          label_id <- which(rownames(explanation) == j)
          labels[label_id] <- paste0(j, " = ", newdata[i, j], 
                                     " (reference category = ", 
                                     levels(object$data[ , j])[1L], ")")
        }
      }
      if (!is.null(rule.col)) {
        text(y = midpoints, x = maxval[1L], labels = labels, cex = cex, 
             col = rule.col, pos = 4L)
      }
    }
  }
  if (intercept) newdata <- cbind(`(Intercept)` = 1L, newdata)
  
  return(list(predictors = newdata[ , req_vars], contribution = t(explanation), predicted.value = preds))
}