


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
#' @param ... Further arguments to be passed to \code{\link{predict.pre}}.
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
#' @references Fokkema, M. (2020). Fitting prediction rule ensembles with R 
#' package pre. \emph{Journal of Statistical Software, 92}(12), 1-30.
#' \doi{10.18637/jss.v092.i12}
#' 
#' Friedman, J. H., & Popescu, B. E. (2008). Predictive learning 
#' via rule ensembles. \emph{The Annals of Applied Statistics, 2}(3), 916-954,
#' \doi{10.1214/07-AOAS148}.
#' @seealso \code{\link{pre}}, \code{\link{interact}} 
#' @export
bsnullinteract <- function(object, nsamp = 10, parallel = FALSE,
                           penalty.par.val = "lambda.1se", verbose = FALSE,
                           ...)
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
                              penalty.par.val = penalty.par.val, ...)
      # step 4: third part of formula 47 of F&P2008:
      # Calculate predictions F_A(x_p):
      F_A_of_x_p <- predict.pre(bs.ens.null, newdata = bsdataset,
                                penalty.par.val = penalty.par.val, ...)
      # step 5: Calculate ytilde of formula 47 of F&P2008:
      ytilde <- F_A_of_x + object$data[bs_inds, object$y_names] - F_A_of_x_p
      # step 6: Build a model using (x,ytilde), using the same procedure as was
      # originally applied to (x,y):
      bsintmodcall$data <- object$data
      bsintmodcall$data[ , all.vars(object$call$formula[[2]])] <- ytilde
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
      F_A_of_x <- predict.pre(bs.ens.null, newdata = object$data, type = "link",
                              ...)
      # step 4: third part of formula 47 of F&P2008:
      # Calculate predictions F_A(x_p):
      F_A_of_x_p <- predict.pre(bs.ens.null, newdata = bsdataset,
                                penalty.par.val = penalty.par.val,
                                type = "link", ...)
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
Hsquaredj <- function(object, varname, k = 10, penalty.par.val = NULL, verbose = FALSE, ...) {
  # Calculate the predicted value F(x) of the full model for each observation:
  preds_x <- predict.pre(object, newdata = object$data, penalty.par.val = penalty.par.val, ...)
  # Calculate the expected value of F_j(x_j), over all observed values x_/j,
  # and the expected value of F_/j(x_/j), over all observed values x_j:
  exp_dataset <- object$data[rep(row.names(object$data),
                                 times = nrow(object$data)),]
  exp_dataset[,varname] <- rep(object$data[ , varname], each = nrow(object$data))
  # using predict.pre for a hughe dataset may lead to errors, so split
  # computations up in k parts:
  exp_dataset$ids <- sample(1:k, nrow(exp_dataset), replace = TRUE)
  for(i in 1:k) {
    if (verbose) cat(".")
    exp_dataset[exp_dataset$ids == i, "yhat"] <- predict.pre(
      object, newdata = exp_dataset[exp_dataset$ids==i,],
      penalty.par.val = penalty.par.val, ...)
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
#' @param gamma Mixing parameter for relaxed fits. See  
#' \code{\link[glmnet]{coef.cv.glmnet}}.
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
#' @param ... Further arguments to be passed to \code{barplot}.
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
#' @references Fokkema, M. (2020). Fitting prediction rule ensembles with R 
#' package pre. \emph{Journal of Statistical Software, 92}(12), 1-30.
#' \doi{10.18637/jss.v092.i12}
#' 
#' Friedman, J. H., & Popescu, B. E. (2008). Predictive learning 
#' via rule ensembles. \emph{The Annals of Applied Statistics, 2}(3), 916-954,
#' \doi{10.1214/07-AOAS148}.
#' @seealso \code{\link{pre}}, \code{\link{bsnullinteract}} 
#' @export
interact <- function(object, varnames = NULL, 
                     penalty.par.val = "lambda.1se", 
                     gamma = NULL, nullmods = NULL,
                     quantprobs = c(.05, .95),
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
      if (is.null(gamma)) {
        Hsquaredj(object = object, varname = varnames[i], k = k,
                  penalty.par.val = penalty.par.val, verbose = verbose)
      } else {
        Hsquaredj(object = object, varname = varnames[i], k = k,
                  penalty.par.val = penalty.par.val, verbose = verbose, gamma = gamma)
      }
    }
    names(H) <- varnames
    if (!is.null(nullmods)) {
      nullH <- foreach::foreach(i = 1:length(varnames), .combine = "cbind", .packages = "pre") %dopar% {
        # Calculate H_j for the bootstrapped null models:
        nullH <- c()
        for(j in 1:length(nullmods)) {
          if (is.null(gamma)) {
            nullH[j] <- Hsquaredj(object = nullmods[[j]], varname = varnames[i],
                                  k = k, penalty.par.val = penalty.par.val,
                                  verbose = verbose)
          } else {
            nullH[j] <- Hsquaredj(object = nullmods[[j]], varname = varnames[i],
                                  k = k, penalty.par.val = penalty.par.val,
                                  verbose = verbose, gamma = gamma)
          }
        }
        nullH
      }
      nullH <- data.frame(nullH, stringsAsFactors = FALSE)
      names(H) <- colnames(nullH) <- varnames
    }
  } else { # if not parallel computation:
    H <- c()
    if (is.null(nullmods)) {
      for (i in 1:length(varnames)) {
        if (is.null(gamma)) {
          H[i] <- Hsquaredj(object = object, varname = varnames[i],
                            k = k, penalty.par.val = penalty.par.val,
                            verbose = verbose)
        } else {
          H[i] <- Hsquaredj(object = object, varname = varnames[i],
                            k = k, penalty.par.val = penalty.par.val,
                            verbose = verbose, gamma = gamma)
        }
      }
    } else { # Compute H values for the training data and bootstrapped null models
      nullH <- data.frame()
      for (i in 1:length(varnames)) {
        if (is.null(gamma)) {
          H[i] <- Hsquaredj(object = object, varname = varnames[i],
                            k = k, penalty.par.val = penalty.par.val,
                            verbose = verbose)
        } else {
          H[i] <- Hsquaredj(object = object, varname = varnames[i],
                            k = k, penalty.par.val = penalty.par.val,
                            verbose = verbose, gamma = gamma)
        }
        for (j in 1:length(nullmods)) {
          if (is.null(gamma)) {
            nullH[j, i] <- Hsquaredj(object = nullmods[[j]], varname = varnames[i],
                                     k = k, penalty.par.val = penalty.par.val,
                                     verbose = verbose)
          } else {
            nullH[j, i] <- Hsquaredj(object = nullmods[[j]], varname = varnames[i],
                                     k = k, penalty.par.val = penalty.par.val,
                                     verbose = verbose, gamma = gamma)
          }
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


