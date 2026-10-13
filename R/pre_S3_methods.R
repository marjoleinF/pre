

#' Print method for objects of class pre
#'
#' \code{print.pre} prints information about the generated prediction rule 
#' ensemble to the command line
#' 
#' @param x An object of class \code{\link{pre}}.
#' @param penalty.par.val character or numeric. Value of the penalty parameter
#' \eqn{\lambda} to be employed for selecting the final ensemble. The default
#' \code{"lambda.1se"} employs the \eqn{\lambda} value within 1 standard
#' error of the minimum cross-validated error. Alternatively, 
#' \code{"lambda.min"} may be specified, to employ the \eqn{\lambda} value
#' with minimum cross-validated error, or a numeric value \eqn{>0} may be 
#' specified, with higher values yielding a sparser ensemble. To evaluate the 
#' trade-off between accuracy and sparsity of the final ensemble, inspect
#' \code{pre_object$glmnet.fit} and \code{plot(pre_object$glmnet.fit)}.
#' @param digits Number of decimal places to print
#' @param ... Further arguments to be passed to 
#' \code{\link[glmnet]{coef.cv.glmnet}}.
#' @return Prints information about the fitted prediction rule ensemble.
#' @details Note that the CV error is estimated with data that was also used 
#' for learning rules and may be too optimistic. Use function \code{\link{cvpre}} to 
#' obtain a more realistic estimate of future prediction error.
#' @examples \donttest{set.seed(42)
#' airq.ens <- pre(Ozone ~ ., data = airquality[complete.cases(airquality),])
#' print(airq.ens)}
#' @method print pre
#' @seealso \code{\link{pre}}, \code{\link{summary.pre}}, \code{\link{plot.pre}}, 
#' \code{\link{coef.pre}}, \code{\link{importance.pre}}, \code{\link{predict.pre}}, 
#' \code{\link{interact}}, \code{\link{cvpre}} 
#' @export
print.pre <- function(x, penalty.par.val = "lambda.1se", 
                      digits = getOption("digits"),
                      ...) {
  
  if (!inherits(x, c("pre", "gpe"))) {
    stop("Argument x should be of class 'pre' (or 'gpe').")
  }
  
  if (!(length(penalty.par.val) == 1L)) {
    stop("Argument penalty.par.val should be a vector of length 1.")
  } else if (!(penalty.par.val == "lambda.min" || 
               penalty.par.val == "lambda.1se" || 
               (is.numeric(penalty.par.val) && penalty.par.val >= 0))) {
    stop("Argument penalty.par.val should be equal to 'lambda.min', 'lambda.1se' or a numeric value >= 0.")
  }
  
  if (!(length(digits) == 1L && digits == as.integer(digits))) {
    stop("Argument digits should be a single integer.")
  }
  
  ## Print summary
  summary.pre(x, penalty.par.val = penalty.par.val, digits = digits, ...)
  
  ## Print coefficients
  coefs <- coef(x, penalty.par.val = penalty.par.val, ...)
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
#' @param ... Further arguments to be passed to \code{\link[glmnet]{coef.cv.glmnet}}.
#' @return Prints information about the fitted prediction rule ensemble.
#' @details Note that the cv error is estimated with data that was also used 
#' for learning rules and may be too optimistic. Use \code{\link{cvpre}} to 
#' obtain a more realistic estimate of future prediction error.
#' @examples \donttest{set.seed(42)
#' airq.ens <- pre(Ozone ~ ., data = airquality[complete.cases(airquality),])
#' summary(airq.ens)}
#' @method summary pre
#' @seealso \code{\link{pre}}, \code{\link{print.pre}}, \code{\link{plot.pre}}, 
#' \code{\link{coef.pre}}, \code{\link{importance.pre}}, \code{\link{predict.pre}}, 
#' \code{\link{interact}}, \code{\link{cvpre}} 
#' @export
summary.pre <- function(object, penalty.par.val = "lambda.1se", 
                        digits = getOption("digits"), ...) {
  
  if (!inherits(object, c("pre", "gpe"))) {
    stop("Argument object should be of class 'pre' (or 'gpe').")
  }
  
  if (!(length(penalty.par.val) == 1L)) {
    stop("Argument penalty.par.val should be a vector of length 1.")
  } else if (!penalty.par.val %in% c("lambda.min", "lambda.1se") && 
             !(is.numeric(penalty.par.val) && penalty.par.val >= 0)) {
    stop("Argument penalty.par.val should be equal to 'lambda.min', 'lambda.1se' or a numeric value >= 0.")
  }
  
  if (!(length(digits) == 1L && digits == as.integer(digits))) {
    stop("Argument digits should be a single integer.")
  }
  
  cl <- match.call()
  
  if (inherits(object$glmnet.fit, "cv.relaxed")) {
    cl$gamma <- eval(cl$gamma)
    
    ## check if gamma value is specified, and if so whether it is a single, proper value
    if (!is.null(cl$gamma)) {
      if (!(length(cl$gamma) == 1L && cl$gamma >= 0 && cl$gamma <= 1)) {
        stop("Argument gamma has been supplied, but should be a single numeric value [0, 1].")
      }
      if (!cl$gamma %in% object$glmnet.fit$relaxed$gamma) {
        stop("Specified gamma value should be one of: ", 
             paste(object$glmnet.fit$relaxed$gamma, sep = ", "),
             ".")
      }
    }
    
    cl$penalty.par.val <- eval(cl$penalty.par.val)
    
    if (penalty.par.val == "lambda.1se") {
      lambda_ind <- object$glmnet.fit$relaxed$index["1se", 1]
      if (is.null(cl$gamma)) {
        gamma_ind <- object$glmnet.fit$relaxed$index["1se", 2L]
      } else {
        gamma_ind <- which(object$glmnet.fit$relaxed$gamma == cl$gamma)
      }
      cat("\nFinal ensemble with cv error within 1se of minimum: \n\n  lambda = ", 
          object$glmnet.fit$relaxed$lambda.1se,
          "\n  gamma = ", object$glmnet.fit$relaxed$gamma[[gamma_ind]])
    }
    if (penalty.par.val == "lambda.min") {
      lambda_ind <- object$glmnet.fit$relaxed$index["min", 1]
      if (is.null(cl$gamma)) {
        gamma_ind <- object$glmnet.fit$relaxed$index["min", 2]
      } else {
        gamma_ind <- which(object$glmnet.fit$relaxed$gamma == cl$gamma)
      }
      cat("Final ensemble with minimum cv error: \n\n  lambda = ", 
          object$glmnet.fit$relaxed$lambda.min,
          "\n  gamma = ", object$glmnet.fit$relaxed$gamma[[gamma_ind]])
    }
    if (is.numeric(penalty.par.val)) {
      if (is.null(cl$gamma)) {
        stop("Relaxed lasso was employed and numeric penalty.par.val was specified; also specify value for gamma.")
      } else {
        if (cl$gamma %in% object$glmnet.fit$relaxed$gamma) {
          gamma_ind <- which(object$glmnet.fit$relaxed$gamma == cl$gamma)
        } else {
          stop("Rule ensemble was fitted with gamma values: ", 
               object$glmnet.fit$relaxed$gamma, 
               ". Argument gamma should specify one of those values.")
        }
        lambda_ind <- which.min(abs(
          object$glmnet.fit$relaxed$statlist[[gamma_ind]]$lambda - penalty.par.val))
        cat("Final ensemble: \n\n  lambda = ", 
            object$glmnet.fit$relaxed$statlist[[gamma_ind]]$lambda[lambda_ind],
            "\n  gamma = ", object$glmnet.fit$relaxed$gamma[gamma_ind])
      }
    }
    cat("\n  number of terms = ", object$glmnet.fit$relaxed$statlist[[gamma_ind]]$nzero[lambda_ind], 
        "\n  mean cv error (se) = ", object$glmnet.fit$relaxed$statlist[[gamma_ind]]$cvm[lambda_ind], 
        " (", object$glmnet.fit$relaxed$statlist[[gamma_ind]]$cvsd[lambda_ind], ")", "\n\n  cv error type : ",
        object$glmnet.fit$name, "\n\n", sep = "")
    
  } else { ## default, non-relaxed lasso was fitted
    
    cl$penalty.par.val <- eval(cl$penalty.par.val)
    
    if (!is.null(cl$gamma)) { 
      warning("A value for gamma was specified, but will be ignored because the rule ensemble was not fit using relax = TRUE")
    }
    
    if (penalty.par.val == "lambda.1se") {
      lambda_ind <- object$glmnet.fit$index["1se", 1L]
      cat("\nFinal ensemble with cv error within 1se of minimum: \n\n  lambda = ", 
          object$glmnet.fit$lambda.1se)
    }
    if (penalty.par.val == "lambda.min") {
      lambda_ind <- object$glmnet.fit$index["min", 1L]
      cat("Final ensemble with minimum cv error: \n\n  lambda = ", 
          object$glmnet.fit$lambda.min)
    }
    if (is.numeric(penalty.par.val)) {
      lambda_ind <- which(abs(object$glmnet.fit$lambda - penalty.par.val) == min(abs(
        object$glmnet.fit$lambda - penalty.par.val)))
      cat("Final ensemble: \n\n  lambda = ", object$glmnet.fit$lambda[lambda_ind])
    }
    cat("\n  number of terms = ", object$glmnet.fit$nzero[lambda_ind], 
        "\n  mean cv error (se) = ", object$glmnet.fit$cvm[lambda_ind], 
        " (", object$glmnet.fit$cvsd[lambda_ind], ")", "\n\n  cv error type : ",
        object$glmnet.fit$name, "\n\n", sep = "")
  }
}









#' Coefficients for the final prediction rule ensemble
#'
#' \code{coef.pre} returns coefficients for prediction rules and linear terms in 
#' the final ensemble
#' 
#' @param object object of class \code{\link{pre}}
#' @inheritParams print.pre
#' @param ... Further arguments to be passed to \code{\link[glmnet]{coef.cv.glmnet}}.
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
#' \code{\link{cvpre}}, \code{\link{importance.pre}}, \code{\link{predict.pre}}, 
#' \code{\link{interact}}, \code{\link{print.pre}} 
#' @export
coef.pre <- function(object, penalty.par.val = "lambda.1se", ...)
{
  
  ## TODO: Add argument on whether learners with zero coefficients should
  ## be included or not
  
  ## check if proper object argument is specified:
  if (!inherits(object, "pre")) {
    stop("Argument object should supply an object of class 'pre'")
  }
  
  ## check if proper penalty.par.val argument is specified:
  if (!(length(penalty.par.val) == 1L)) {
    stop("Argument penalty.par.val should be a vector of length 1.")
  } else if (!(penalty.par.val == "lambda.min" || 
               penalty.par.val == "lambda.1se" || 
               (is.numeric(penalty.par.val) && penalty.par.val >= 0))) {
    stop("Argument penalty.par.val should be equal to 'lambda.min', 'lambda.1se' or a numeric value >= 0")
  }
  
  ## check if gamma value is specified, and if so whether it is a single value
  gamma <- eval.parent(match.call()[["gamma"]])
  if (!is.null(gamma)) {  
    if (is.null(object$glmnet.fit$relaxed)) {
      warning("A gamma value was specified, but pre object was not fitted using relax = TRUE. Specified gamma will be ignored or an error may occur.")
    } else {
      #gamma <- eval.parent(parse(text = cl$gamma))
      if (!(length(gamma) == 1L && gamma >= 0 && gamma <= 1)) {
        stop("Argument gamma has been supplied, but should be a single numeric value [0, 1].")
      }
      if (!gamma %in% object$glmnet.fit$relaxed$gamma) {
        stop("Specified gamma value should be one of ", object$glmnet.fit$relaxed$gamma)
      }
      if (is.null(object$glmnet.fit$relaxed)) {
        warning("A gamma value was specified, but the pre ensemble was not fitted using relax = TRUE. The gamma value will be ignored.")
      }
    }
  }
  
  if (object$family %in% c("gaussian", "binomial", "poisson", "cox")) {
    coefs <- as(coef(object$glmnet.fit, s = penalty.par.val, ...), 
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
#' @param newdata optional \code{data.frame} of new (test) observations, including all
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
#' predict(airq.ens, newdata = airquality[complete.cases(airquality),][-train,])
#' 
#' }
#' @import Matrix
#' @method predict pre
#' @seealso \code{\link{pre}}, \code{\link{plot.pre}}, 
#' \code{\link{coef.pre}}, \code{\link{importance.pre}}, \code{\link{cvpre}}, 
#' \code{\link{interact}}, \code{\link{print.pre}}, 
#' \code{\link[glmnet]{predict.cv.glmnet}}
#' @export
predict.pre <- function(object, newdata = NULL, type = "link",
                        penalty.par.val = "lambda.1se", ...)
{
  ## Check if proper object argument is specified:
  if (!inherits(object, "pre")) {
    stop("Argument object should supply an object of class 'pre'")
  }
  
  ## check if proper type argument is specified:
  if (length(type) != 1L || !is.character(type)) {
    stop("Argument type should be a character vector of length 1")
  }
  
  ## check if proper penalty.par.val argument is specified:
  penalty.par.val <- eval(penalty.par.val)  
  if (!(length(penalty.par.val) == 1L)) {
    stop("Argument penalty.par.val should be a vector of length 1.")
  } else if (!(penalty.par.val%in% c("lambda.min", "lambda.1se")) && 
             !(is.numeric(penalty.par.val) && penalty.par.val >= 0)) {
    stop("Argument penalty.par.val should be equal to 'lambda.min', 'lambda.1se' or a numeric value >= 0")
  }
  
  ## check if gamma value is specified, and if so whether it is a single value
  cl <- match.call()
  if (!is.null(cl$gamma)) {
    cl$gamma <- eval(cl$gamma)
    if (!(length(cl$gamma) == 1L && cl$gamma >= 0 && cl$gamma <= 1)) {
      stop("Argument gamma has been supplied, but should be a single numeric value [0, 1].")
    }
    if (!cl$gamma %in% object$glmnet.fit$relaxed$gamma) {
      stop("Specified gamma value should be one of ", object$glmnet.fit$relaxed$gamma)
    }
    if (is.null(object$glmnet.fit$relaxed)) {
      warning("A gamma value was specified, but the pre ensemble was not fit using relax = TRUE. The gamma value will be ignored.")
    }
  }
  
  ## Construct data matrix for prediction
  if (is.null(newdata)) {
    
    newdata <- object$modmat
    
  } else {
    
    if (inherits(newdata, c("tbl_df", "tbl"))) newdata <- as.data.frame(newdata)
    
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
      if (any(ordered_names <- sapply(newdata, is.ordered))) {
        ordered_names <- names(newdata)[ordered_names]
        newdata[ , ordered_names] <- sapply(newdata[ , ordered_names], as.numeric)
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
      confirmatory = eval(object$call$confirmatory))
    
    newdata <- newdata$x
  }
  
  ## Get predictions:
  if (object$family %in% c("gaussian", "binomial", "poisson", "cox")) {
    
    preds <- predict(object$glmnet.fit, newx = newdata, 
                     s = penalty.par.val, type = type, ...)[ , 1L]
  } else if (object$family %in% c("mgaussian", "multinomial")) {
    if (object$family == "multinomial" && type == "class") {
      preds <- predict(object$glmnet.fit, newx = newdata, 
                       s = penalty.par.val, type = type, ...)[ , 1L]
    } else {
      preds <- predict(object$glmnet.fit, newx = newdata, 
                       s = penalty.par.val, type = type, ...)[ , , 1L]
    }
  }
  return(preds)
}






##' @export
importance <- function(x, ...)  UseMethod("importance")


#' Calculate importances of baselearners and input variables in a prediction 
#' rule ensemble (pre)
#'
#' \code{importance.pre} calculates importances for rules, linear terms and input
#' variables in the prediction rule ensemble (pre), and creates a bar plot 
#' of variable importances.
#'
#' @param x an object of class \code{\link{pre}}
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
#' @param gamma Mixing parameter for relaxed fits. See  
#' \code{\link[glmnet]{coef.cv.glmnet}}.
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
#' @references Fokkema, M. (2020). Fitting prediction rule ensembles with R 
#' package pre. \emph{Journal of Statistical Software, 92}(12), 1-30.
#' \doi{10.18637/jss.v092.i12}
#' 
#' Fokkema, M. & Strobl, C. (2020). Fitting prediction rule ensembles to psychological 
#' research data: An introduction and tutorial. \emph{Psychological Methods 25}(5), 
#' 636-652. \doi{10.1037/met0000256}, \url{https://arxiv.org/abs/1907.05302}
#' 
#' Friedman, J. H., & Popescu, B. E. (2008). Predictive learning 
#' via rule ensembles. \emph{The Annals of Applied Statistics, 2}(3), 916-954 
#' \doi{10.1214/07-AOAS148}.
#' @seealso \code{\link{pre}}
#' @export
#' @method importance pre
#' @aliases importance
importance.pre <- function(x, standardize = FALSE, global = TRUE,
                           penalty.par.val = "lambda.1se", gamma = NULL,
                           quantprobs = c(.75, 1),
                           round = NA, plot = TRUE, ylab = "Importance",
                           main = "Variable importances", abbreviate = 10L,
                           diag.xlab = TRUE, diag.xlab.hor = 0, diag.xlab.vert = 2,
                           cex.axis = 1, legend = "topright", ...)
{
  
  if (!inherits(x, what = "pre")) {
    stop("Specified object is not of class 'pre'.")
  }
  
  if (!global) {
    if (x$family %in% c("mgaussian", "multinomial")) {
      warning("Local importances cannot be calculated for multivariate and multinomial outcomes. Global importances will be returned.")
      global <- TRUE 
    }
  }
  
  if (standardize && x$family %in% c("multinomial", "binomial", "cox")) {
    warning("Standardized importances cannot be calculated for binary, multinomial or survival responses. Unstandardized importances will be returned.")
    standardize <- FALSE
  }
  
  ## Step 1: Calculate the importances of the base learners:
  
  ## get coefficients:
  if (is.null(gamma)) {
    coefs <- coef(x, penalty.par.val = penalty.par.val)
  } else {
    coefs <- coef(x, penalty.par.val = penalty.par.val, gamma = gamma)    
  }
  if (x$family %in% c("mgaussian", "multinomial")) {
    coef_inds <- names(coefs)[!names(coefs) %in% c("rule", "description")]
  }
  
  ## continue only when there are nonzero terms besides intercept:
  if ((x$family %in% c("gaussian", "binomial", "poisson") && 
       sum(coefs$coefficient != 0) > 1L ) || 
      (x$family %in% c("mgaussian", "multinomial") && 
       sum(rowSums(coefs[,coef_inds]) != 0) > 1L) ||
      (x$family == "cox" && sum(coefs$coefficient != 0) > 0)) { 
    ## give factors a description:
    if (any(is.na(coefs$description))) {
      coefs$description[is.na(coefs$description)] <-
        paste0(as.character(coefs$rule)[is.na(coefs$description)], " ")
    }
    coefs <- coefs[order(coefs$rule), ]
    
    ## Get SDs for every baselearner:
    if (global) {
      ## Get SDs (x$x_scales should be used to get correct SDs for linear terms)
      if (x$family == "cox") {
        sds <- apply(x$modmat, 2, sd, na.rm = TRUE)  
      } else {
        sds <- c(0, apply(x$modmat, 2, sd, na.rm = TRUE))          
      }
      if (standardize) {
        if (x$family == "mgaussian") {
          sd_y <- sapply(x$data[ , x$y_names], sd)
        } else if (x$family %in% c("gaussian", "poisson")) { 
          sd_y <- sd(as.numeric(x$data[ , x$y_names]))
        }
      }
    } else {
      preds <- predict.pre(x, newdata = x$data, type = "response",
                           penalty.par.val = penalty.par.val, ...)
      local_modmat <- x$modmat[preds >= quantile(preds, probs = quantprobs[1]) &
                                 preds <= quantile(preds, probs = quantprobs[2]),]
      if (nrow(local_modmat) < 2) {stop("Selected subregion contains less than 2 observations, importances cannot be calculated")}
      ## x$x_scales should be used to get correct SDs for linear terms:
      if (x$family == "cox") {
        ## cox prop hazard model has no intercept, so should be omitted
        sds <- apply(local_modmat, 2, sd, na.rm = TRUE)
      } else {
        sds <- c(0, apply(local_modmat, 2, sd, na.rm = TRUE))
      }
      
      if (standardize) {
        sd_y <- sd(x$data[preds >= quantile(preds, probs = quantprobs[1]) & 
                            preds <= quantile(preds, probs = quantprobs[2]),
                          x$y_names])
      }
    }
    
    ## Check if there are any " ` " marks in sd names, if so remove:
    if (any(grepl("`", names(sds), fixed = TRUE))) {
      names(sds) <- gsub("`", "", names(sds), fixed = TRUE)
    }
    
    if(x$normalize) {
      sds[names(x$x_scales)] <- sds[names(x$x_scales)] * x$x_scales
    }
    
    if (x$family != "cox") {
      names(sds)[1] <- "(Intercept)"
    }
    
    sds <- sds[order(names(sds))]
    
    if (any(names(sds) != coefs$rule)) {
      warning("There seems to be a problem with the ordering or size of the coefficient and sd vectors. Importances cannot be calculated.")
    }
    
    ## baselearner importance is given by abs(coef*SD) (F&P section 6):
    if (x$family %in% c("multinomial", "mgaussian")) {
      baseimps <- data.frame(coefs, sd = sds)
      baseimps[,gsub("coefficient", "importance", coef_inds)] <- abs(sapply(baseimps[,coef_inds], function(x) x*sds))
    } else {
      baseimps <- data.frame(coefs, sd = sds, imp = abs(coefs$coefficient)*sds)
    }
    
    if (standardize) {
      if (x$family == "mgaussian") {
        for (i in gsub("coefficient", "importance", coef_inds)) {
          baseimps[,i] <- baseimps[,i] / sd_y[gsub("importance.", "", i)]
        }
      } else if (x$family %in% c("gaussian", "poisson")) {
        baseimps$imp <- baseimps$imp / sd_y
      }
    }
    
    ## Remove nonzero terms:
    if (x$family %in% c("mgaussian", "multinomial")) {
      baseimps <- baseimps[rowSums(baseimps[,coef_inds]) != 0, ]   
    } else {
      baseimps <- baseimps[baseimps$coefficient != 0,]
    }
    
    ## Omit intercept:
    baseimps <- baseimps[baseimps$description != "1",]
    
    ## Calculate the number of conditions in each rule:
    baseimps$nterms <- NA
    for(i in 1:nrow(baseimps)) {
      ## If there is " & " in description, there are at least 2 conditions/variables 
      ## in the base learner:
      if (grepl(" & ", baseimps$description[i])) {
        baseimps$nterms[i] <- length(gregexpr("&", baseimps$description)[[i]]) + 1L
      } else {
        baseimps$nterms[i] <- 1L # if not, the number of terms equals 1
      }
    }
    
    ## if no winsorizing is performed, descriptions look different then with winsorizing
    ##    so add a temporary space AFTER description
    if (!is.null(x$call$winsfrac)) {
      if (x$call$winsfrac == 0) {
        linear_term_ids <- which(baseimps$rule == baseimps$description)
        for (i in linear_term_ids) {
          baseimps$description[i] <- paste0(baseimps$description[i], " ")
        }
      }
    }
    
    
    ## Step 2: Calculate variable importances:
    if (x$family %in% c("mgaussian", "multinomial")) {
      varimps <- data.frame(varname = x$x_names, stringsAsFactors = FALSE)
      varimps[,gsub("coefficient", "importance", coef_inds)] <- 0
    } else {
      varimps <- data.frame(varname = x$x_names, imp = 0,
                            stringsAsFactors = FALSE)
    }
    
    for(i in 1:nrow(varimps)) {
      ## Get imps from rules and linear functions
      for(j in 1:nrow(baseimps)) {
        
        ## if the variable name appears in the description (of rule or linear term):
        ##   (Note: EXACT matches are needed, so 1) there should be a space before 
        ##     and after the variable name in the rule and thus 2) there should be 
        ##     a space added before the description of the rule)
        if (grepl(paste0(" ", varimps$varname[i], " "), paste0(" ", baseimps$description[j]))) {
          ## Count the number of times it appears in the rule
          n_occ <- length(gregexpr(paste0(" ", varimps$varname[i], " "),
                                   paste0(" ", baseimps$description[j]), fixed = TRUE)[[1]])
          ## Add to the importance of the variable
          if (x$family %in% c("mgaussian", "multinomial")) {
            varimps[i, gsub("coefficient", "importance", coef_inds)] <- 
              varimps[i, gsub("coefficient", "importance", coef_inds)] + 
              (n_occ * baseimps[j, gsub("coefficient", "importance", coef_inds)] / baseimps$nterms[j])
          } else {
            varimps$imp[i] <- varimps$imp[i] + (n_occ * baseimps$imp[j] / baseimps$nterms[j])
          }
        }
      }
      
      ## Get imps for factors
      if (is.factor(x$data[ , varimps$varname[i]])) { # check if variable is a factor and add importance
        # !is.ordered(x$data[ , varimps$varname[i]])) {
        ## Sum baseimps$imp for which baseimps$rule has varimps$varname[i] as _part_ of its name
        if (x$family %in% c("mgaussian", "multinomial")) {
          varimps[i, gsub("coefficient", "importance", coef_inds)] <-
            varimps[i, gsub("coefficient", "importance", coef_inds)] +
            colSums(baseimps[grepl(varimps$varname[i], baseimps$rule, fixed = TRUE), 
                             gsub("coefficient", "importance", coef_inds)])
        } else { # otherwise, simply add importance(s)
          varimps$imp[i] <- varimps$imp[i] + 
            sum(baseimps$imp[grepl(varimps$varname[i], baseimps$rule, fixed = TRUE)])
        }
      } 
    }
    
    
    ## Step 3: Return (and plot) importances:
    
    if (x$family %in% c("mgaussian", "multinomial")) {
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
      if (x$family %in% c("mgaussian", "multinomial")) {
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
                  cex.axis = cex.axis, ...)
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
                  main = main, cex.axis = cex.axis, ...)
        }
      }
    }
    
    if (!is.na(round)) {
      baseimps[,sapply(baseimps, is.numeric)] <- round(baseimps[,sapply(baseimps, is.numeric)], digits = round)
      varimps[,sapply(varimps, is.numeric)] <- round(varimps[,sapply(varimps, is.numeric)], digits = round)
    }
    
    if (x$family %in% c("mgaussian","multinomial")) {
      keep <- c("rule", "description", gsub("coefficient", "importance", coef_inds), 
                coef_inds, "sd")
    } else {
      keep <- c("rule", "description", "imp", "coefficient", "sd")
    }
    
    baseimps <- data.frame(baseimps[, keep], stringsAsFactors = FALSE)
    row.names(baseimps) <- row.names(varimps) <- NULL
    ## Remove added space AFTER description if winsorizing was performed
    if (!is.null(x$call$winsfrac)) {
      if (x$call$winsfrac == 0L) {
        for (i in linear_term_ids) {
          baseimps$description[i] <- substring(baseimps$description[i], first = 2L)
        }
      }
    }
    
    return(invisible(list(varimps = varimps, baseimps = baseimps)))
    
  } else {
    warning("No non-zero terms in the ensemble. All importances are zero.")
    return(invisible(NULL))
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
#' \code{\link{importance.pre}}.
#' @param gamma Mixing parameter for relaxed fits. See  
#' \code{\link[glmnet]{coef.cv.glmnet}}.
#' @param ... Arguments to be passed to \code{\link[grid]{gpar}}.
#' @examples
#' \donttest{set.seed(42)
#'  airq.ens <- pre(Ozone ~ ., data = airquality[complete.cases(airquality),])
#'  plot(airq.ens)}
#' @seealso \code{\link{pre}}, \code{\link{print.pre}}
#' @method plot pre
#' @export
plot.pre <- function(x, penalty.par.val = "lambda.1se", gamma = NULL,
                     linear.terms = TRUE, 
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
    if (is.null(gamma)) {
      coefs <- coef(x, penalty.par.val = penalty.par.val)
    } else {
      coefs <- coef(x, penalty.par.val = penalty.par.val, gamma = gamma)
    }
    nonzeroterms <- coefs[rowSums(coefs[,!names(coefs) %in% c("rule", "description")]) != 0,]
    if ("(Intercept)" %in% nonzeroterms$rule) {
      intercept <- nonzeroterms[which(nonzeroterms$rule == "(Intercept)"), "coefficient"] # may be needed for plotting linear terms later      
      nonzeroterms <- nonzeroterms[-which(nonzeroterms$rule == "(Intercept)"), ] # omit intercept
    }
  } else {
    if (is.null(gamma)) {
      nonzeroterms <- importance(x, plot = FALSE, global = TRUE, 
                                 penalty.par.val = penalty.par.val, 
                                 standardize = standardize)$baseimps
    } else {
      nonzeroterms <- importance(x, plot = FALSE, global = TRUE, 
                                 penalty.par.val = penalty.par.val, 
                                 standardize = standardize, gamma = gamma)$baseimps
    }
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
          levels(treeplotdata[ , j]) <- c(
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
          ht <- nlevels(x$data [ ,x$y_names])
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