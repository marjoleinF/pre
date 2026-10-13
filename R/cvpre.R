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
#' @param ... Further arguments to be passed to \code{\link{predict.pre}}.
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
#' \code{\link{coef.pre}}, \code{\link{importance.pre}}, \code{\link{predict.pre}}, 
#' \code{\link{interact}}, \code{\link{print.pre}} 
#' @export
cvpre <- function(object, k = 10, penalty.par.val = "lambda.1se", pclass = .5, 
                  foldids = NULL, verbose = FALSE, parallel = FALSE,
                  print = TRUE, ...) {
  
  ## check if proper object argument is specified:
  if (!inherits(object, "pre")) {
    stop("Argument object should supply an object of class 'pre'")
  }
  
  ## Check if proper k argument is specified:  
  if (!(length(k) == 1L && k == as.integer(k))) {
    stop("Argument k should be a single positive integer.")
  }  
  
  ## Check if proper verbose argument is specified:  
  if (!(is.logical(verbose) && length(verbose) == 1L)) {
    stop("Argument verbose should be TRUE or FALSE.")
  }  
  
  ## check if pclass is a numeric vector of length 1, and <= 1 and > 0
  if (!(is.numeric(pclass) && length(pclass) == 1L && pclass <= 1 && pclass > 0)) {
    stop("Argument verbose should be TRUE or FALSE.")
  }  
  
  ## check if proper penalty.par.val argument is specified:
  if (!(length(penalty.par.val) == 1L)) {
    stop("Argument penalty.par.val should be a numeric vector of length 1.")
  } else if (!(penalty.par.val == "lambda.min" || 
               penalty.par.val == "lambda.1se" || 
               (is.numeric(penalty.par.val) && penalty.par.val >= 0))) {
    stop("Argument penalty.par.val should be equal to 'lambda.min', 'lambda.1se' or a numeric value >= 0")
  }
  
  ## check if proper parallel argument is specified:
  if (!(is.logical(parallel) && length(parallel) == 1L)) {
    stop("Argument parallel should be TRUE or FALSE")
  }
  
  ## check if proper foldids argument is specified:
  if (!is.null(foldids)) {
    if (length(foldids) != nrow(object$data)) {
      stop("Argument foldids has length ", length(foldids), ", but should have length ", nrow(object$data), ".")
    } else if (!all.equal(foldids, as.integer(foldids))) {
      stop("Argument foldids should be an integer vector, but is not.")
    } 
  }
  
  ## Set up fold-ids, seeds and object for collecting CV predictions:
  if (is.null(foldids)) {
    foldids <- sample(rep(1:k, length.out = nrow(object$data)), 
                      size = nrow(object$data), replace = FALSE)
  }
  seeds <- sample(k*99, size = k)
  y_ncol <- ifelse(object$family == "multinomial", 
                   nlevels(object$data[ , object$y_names]), 
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
              penalty.par.val = penalty.par.val, ...)
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
        type = "response", penalty.par.val = penalty.par.val, ...)
      
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