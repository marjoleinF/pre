#' Get the optimal lambda and gamma parameter values for an ensemble of given size
#'
#' Function \code{prune_pre} returns the optimal values of lambda and gamma for
#' the requested ensemble size.
#' 
#' @param object an object of class \code{\link{pre}} that was fit using the relaxed lasso.
#' If an object of class \code{\link{pre}} is specified that was not fit using the relaxed lasso,
#' an error will be printed.
#' @param nonzero maximum number of terms to retain.
#' @param plusminus number of terms above and below \code{nonzero} for which CV results will be printed.
#' 
#' @return The lambda and gamma values that yield optimal predictive accuracy for the specified
#' number of terms. These are invisibly returned, see Examples on how to use them. A sentence
#' describing what the optimal values are is printed to the command line, with an overview of
#' the performance (in terms of cross-validated accuracy and the number of terms retained) of 
#' lambda values near the optimum. If the specified number of terms to retain is lower than
#' what would be obtained using the \code{lambda.min} or \code{lambda.1se} criterion, a warning
#' will also be printed.
#' @examples \donttest{
#' ## Fit a rule ensemble to predict Ozone concentration
#' airq <- airquality[complete.cases(airquality), ]
#' set.seed(42)
#' airq.ens <- pre(Ozone ~ ., data = airq, relax = TRUE)
#' 
#' ## Inspect the result (default lambda.1se criterion)
#' airq.ens
#' 
#' ## Inspect the lambda path 
#' ## (lower x-axis gives lambda values, upper x-axis corresponding no. of non-zero terms)
#' \dontrun{plot(airq.ens$glmnet.fit)}
#' 
#' ## Accuracy still quite good with only 5 terms, obtain corresponding parameter values
#' opt_pars <- prune_pre(airq.ens, nonzero = 5)
#' opt_pars
#' 
#' ## Use the parameter values for interpretation and prediction, e.g.
#' predict(airq.ens, newdat = airq[c(22, 33), ], penalty = opt_pars$lambda, gamma = opt_pars$gamma)
#' summary(airq.ens, penalty = opt_pars$lambda, gamma = opt_pars$gamma)
#' print(airq.ens, penalty = opt_pars$lambda, gamma = opt_pars$gamma)
#' }
#' @seealso \code{\link{pre}}
#' @export
prune_pre <- function(object, nonzero, plusminus = 3) {
  
  ## Check if nonzero occurs in lambda path, otherwise warn and take nearest number
  if (!nonzero %in% object$glmnet.fit$nzero) {
    nonzero <- object$glmnet.fit$nzero[names(which.min((nonzero - object$glmnet.fit$nzero)^2))[1L]]
    warning(paste0("Specified value of nonzero argument (", nonzero, ") does not occur on the lambda path. Results are returned for", nonzero, "non-zero terms, instead."), 
            immediate. = TRUE)
  }
  
  ## For relaxed lasso
  if (!is.null(object$call$relax) && object$call$relax) {
    
    ## warn if non-optimal sparsity and accuracy are requested
    if (object$glmnet.fit$relaxed$nzero.1se < nonzero) {
      warning("The requested number of non-zero terms (", nonzero, ") is larger than the ", object$glmnet.fit$relaxed$nzero.1se, 
              " non-zero terms retained with the lambda.1se criterion. Both complexity and cross-validated error will likely be lower (better) when using the lambda.1se criterion!", 
              immediate. = TRUE)
    }
    if (object$glmnet.fit$relaxed$nzero.min < nonzero) {
      warning("The requested number of non-zero terms (", nonzero, ") is larger than the ", object$glmnet.fit$relaxed$nzero.min, 
              " non-zero terms retained with the lambda.min criterion. Both complexity and cross-validated error will likely be lower (better) using the default lambda.min criterion!", 
              immediate. = TRUE)
    } 
    ## Get optimal lambda and gamma values for requested number of nonzero terms
    gammas <- names(object$glmnet.fit$relaxed$statlist)
    optimal_gammas <- gammas[apply(sapply(gammas, function(x) object$glmnet.fit$relaxed$statlist[[x]]$cvm), 1L, which.min)]
    df <- data.frame(lambda = object$glmnet.fit$lambda,
                     number_of_nonzero_terms = object$glmnet.fit$nzero,
                     optimal_gamma = optimal_gammas,
                     mean_cv_error = NA)
    svals <- as.numeric(substring(rownames(df)[df$number_of_nonzero_terms == nonzero], first = 2))
    svals <- (min(svals) - plusminus):(max(svals) + plusminus)
    svals <- paste0("s", svals[svals %in% 0:(nrow(df)+1)])
    df <- df[svals, ]
    for (i in rownames(df)) {
      df[i, "mean_cv_error"] <- data.frame(object$glmnet.fit$relaxed$statlist[[df[i, "optimal_gamma"]]])[i, "cvm"]
    }
    df$optimal_gamma <- substring(df$optimal_gamma, first = 3)
    min_cvm <- which.min(df$mean_cv_error[df$number_of_nonzero_terms == nonzero])
    lambda <- df$lambda[df$number_of_nonzero_terms == nonzero][min_cvm]
    gamma <- as.numeric(df$optimal_gamma[df$number_of_nonzero_terms == nonzero][min_cvm])
    cat(paste0("The best ensemble with ", nonzero, " non-zero terms is obtained with a lambda value of ", 
               round(lambda, digits= 6), " and a gamma value of ", gamma, ".\n\n"))
    cat(paste0("Overview of performance of ensembles selected with the nearest lambda values:\n"))
    print(df)
    invisible(list(lambda = lambda, gamma = gamma))
  } else { ## For non-relaxed lasso
    stop("For obtaining an ensemble with a pre-specified number of non-zero terms, use of the relaxed lasso is strongly recommended. The specified ensemble was fit using standard lasso. Please refit the original ensemble using function pre() and additionally specifying relax = TRUE.")
  }
}

