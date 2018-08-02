#' Print a General Prediction Ensemble (gpe)
#' @export
#' @method print gpe
#' @inheritParams print.pre
#' @param x An object of class \code{\link{gpe}}.
#' @seealso \code{\link{gpe}} \code{\link{print.pre}}
print.gpe <- function(
  x, penalty.par.val = "lambda.1se", digits = getOption("digits"), ...){
  out <- print.pre(x, penalty.par.val, digits, ...)
  
  # Make comment about abbreviation of hinge function
  if(any(grepl("^eTerm", out$description)))
    cat("\n  'h' in the 'eTerm' indicates the hinge function\n")
  
  invisible(out)
}


#' Summary method for a General Prediction Ensemble (gpe)
#'
#' \code{summary.gpe} prints information about the generated ensemble 
#' to the command line
#' 
#' @param object An object of class \code{\link{gpe}}.
#' @inheritParams print.pre
#' @return Prints information about the fitted ensemble.
#' @details Note that the cv error is estimated with data that was also used 
#' for learning rules and may be too optimistic.
#' @export
#' @method summary gpe
#' @seealso \code{\link{gpe}}, \code{\link{print.gpe}}, 
#' \code{\link{coef.gpe}}, \code{\link{predict.gpe}}
summary.gpe <- function(object, penalty.par.val = "lambda.1se", ...) {
  
  if (class(object) != "gpe") {
    stop("Argument 'object' should be of class 'gpe'.")
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




#' @title Coefficients for a General Prediction Ensemble (gpe)
#' 
#' @description coef function for \code{\link{gpe}}
#' @export
#' @method coef gpe
#' @inheritParams coef.pre
#' @seealso \code{\link{coef.pre}}
coef.gpe <- function(object, penalty.par.val = "lambda.1se", ...)
{
  coefs <- as(coef.glmnet(object$glmnet.fit, s = penalty.par.val, ...), 
              Class = "matrix")
  colnames(coefs) <- "coefficient"
  
  # Change the row.names
  rownames(coefs) <- gpe_pretty_labels(row.names(coefs))
  
  # Return
  coefs <- coefs[coefs != 0, , drop = FALSE]
  data.frame(description = row.names(coefs), coefs)
}

gpe_pretty_labels <- function(term_labels){
  out <- term_labels
  
  # Handle rTerm
  is_rTerm <- grepl("^rTerm\\(", out)
  out[is_rTerm] <- str_replace_all(out[is_rTerm], "(^rTerm\\()|(\\)$)", "")
  
  # Handle lTerm. We don't do anything due to as we want the scales etc
  # is_lTerm <- grepl("^lTerm\\(", out)
  # out[is_lTerm] <- str_replace_all(out[is_lTerm], "(^lTerm\\()|(\\)$)", "")
  
  # Handle eTerm
  is_eTerm <- grepl("^eTerm\\(", out)
  # Replace match group 2 with 'h(' and 4 with ')'
  expr <- "(eTerm\\(.*)(pmax\\()(.*)(\\,\\ 0\\))(.*)"
  expr_replace <- "\\1h(\\3)\\5"
  while(any(grepl(expr, out[is_eTerm])))
    out[is_eTerm] <- gsub(expr, expr_replace, out[is_eTerm], perl = TRUE)
  
  out
}

#' @title Predicted values based on gpe ensemble
#' 
#' @description 
#' Predict function for \code{\link{gpe}}
#' 
#' @param object of class \code{\link{gpe}}
#' @param newdata optional new data to compute predictions for
#' @param type argument passed to \code{\link{predict.cv.glmnet}}
#' @param penalty.par.val argument passed to \code{s} argument of \code{\link{predict.cv.glmnet}}
#' @param ... Unused
#' 
#' @details 
#' The initial training data is used if \code{newdata = NULL}.
#' 
#' @seealso 
#' \code{\link{gpe}}
#' 
#' @export
#' @method predict gpe
predict.gpe <- function(
  object, newdata = NULL, type = "link",
  penalty.par.val = "lambda.1se", ...){
  if(is.null(newdata)){
    if(is.null(object$data))
      stop("Predict called with no new object and no saved data with gpe")
    newdata <- object$data
  }
  
  X <- model.matrix(object$terms, newdata)
  predict(object$glmnet.fit, newx = X, s = penalty.par.val, type = type)
}