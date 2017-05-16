#' @rdname print.pre
#' @export
#' @method print gpe
print.gpe <- function(x, penalty.par.val = "lambda.1se", ...){
  print.pre(x, penalty.par.val, ...)
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