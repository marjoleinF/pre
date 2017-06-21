#' @rdname print.pre
#' @export
#' @method print gpe
print.gpe <- function(
  x, penalty.par.val = "lambda.1se", digits = getOption("digits"), ...){
  out <- print.pre(x, penalty.par.val, digits, ...)
  
  # Make comment about abbreviation of hinge function
  if(any(grepl("^eTerm", out$description)))
    cat("\n  'h' in the 'eTerm' indicates the hinge function\n")
  
  invisible(out)
}

#' @rdname coef.pre
#' @export
#' @method coef gpe
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

#' @title Predicted Values Based on gpe Ensemble
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