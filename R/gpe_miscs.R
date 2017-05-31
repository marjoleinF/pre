#' @rdname print.pre
#' @export
#' @method print gpe
print.gpe <- function(
  x, penalty.par.val = "lambda.1se", digits = getOption("digits"), ...){
  print.pre(x, penalty.par.val, digits, ...)
}

#' @rdname coef.pre
#' @export
#' @method coef gpe
coef.gpe <- function(object, penalty.par.val = "lambda.1se", ...)
{
  coefs <- as(coef.glmnet(object$glmnet.fit, s = penalty.par.val, ...), 
              Class = "matrix")
  colnames(coefs) <- "coefficient"
  coefs <- coefs[coefs != 0, ,drop = FALSE]
  data.frame(description = row.names(coefs), coefs)
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