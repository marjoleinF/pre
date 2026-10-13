#' Create partial dependence plot for a single variable in a prediction rule 
#' ensemble (pre)
#'
#' \code{singleplot} creates a partial dependence plot, which shows the effect of
#' a predictor variable on the ensemble's predictions. Note that plotting partial 
#' dependence is computationally intensive. Computation time will increase fast 
#' with increasing numbers of observations and variables. For large 
#' datasets, package `plotmo` (Milborrow, 2019) provides more efficient functions 
#' for plotting partial dependence and also supports `pre` models. 
#'
#' @param object an object of class \code{\link{pre}}.
#' @param gamma Mixing parameter for relaxed fits. See  
#' \code{\link[glmnet]{coef.cv.glmnet}}.
#' @param varname character vector of length one, specifying the variable for
#' which the partial dependence plot should be created. Note that \code{varname}
#' should correspond to the variable as described in the model formula used
#' to generate the ensemble (i.e., including functions applied to the variable).
#' @inheritParams print.pre
#' @param response numeric vector of length 1. Only relevant for multivariate gaussian 
#' and multinomial responses. If \code{NULL} (default), PDPs for all response 
#' variables or categories will be produced. A single integer can be specified, 
#' indicating for which response variable or category PDPs should be produced.  
#' @param nvals optional numeric vector of length one. For how many values of x
#' should the partial dependence plot be created?
#' @param type character string. Type of prediction to be plotted on y-axis.
#' \code{type = "response"} gives fitted values for continuous outputs and
#' fitted probabilities for nominal outputs. \code{type = "link"} gives fitted
#' values for continuous outputs and linear predictor values for nominal outputs.
#' @param ylab character. Label to be printed on the y-axis, defaults to the response
#' variable name(s).
#' @param xlab character. Label to be printed on the x-axis. If \code{NULL},
#' the supplied \code{varname} will be printed on the x-axis.
#' @param newdata Optional \code{data.frame} in which to look for variables 
#' with which to predict. If \code{NULL} (the default), the \code{data.frame} used to fit the 
#' original ensemble will be used. Smaller subsets of the original data can
#' be specified to (substantially) reduce computation time. See Details.
#' @param rug logical. Should a rug be plotted on the x-axis, representing the 
#' location of observed datapoint? Note that the rug will only show where values
#' have been observed, not their frequency/density. 
#' @param ... Further arguments to be passed to 
#' \code{\link[graphics]{plot.default}}.
#' @return A 1D partial dependence plot will be plotted. Invisibly, partial dependence
#' values are returned.  
#' @details By default, a partial dependence plot will be created for each unique
#' observed value of the specified predictor variable. See also section 8.1 of 
#' Friedman & Popescu (2008).
#' 
#' When the number of unique observed values is large, partial dependence functions
#' can take a very long time to compute. Specifying the \code{nvals} argument 
#' can substantially reduce computation time. When the
#' \code{nvals} argument is supplied, values for the minimum, maximum, and \code{(nvals - 2)}
#' intermediate values of the predictor variable will be plotted. Note that \code{nvals}
#' can be specified only for numeric and ordered input variables. If the plot is
#' requested for a nominal input variable, the \code{nvals} argument will be
#' ignored and a warning printed.
#' 
#' Alternatively, \code{newdata} can be specified to provide a different (smaller) 
#' set of observations to compute partial dependence over.
#' If \code{mi_pre} was used to derive the original rule ensemble, 
#' function \code{mean_mi} can be used for this.
#' 
#' @references Friedman, J. H., & Popescu, B. E. (2008). Predictive learning 
#' via rule ensembles. \emph{The Annals of Applied Statistics, 2}(3), 916-954.
#' 
#' Milborrow, S. (2019). plotmo: Plot a model's residuals, response, and partial 
#' dependence plots. \url{https://CRAN.R-project.org/package=plotmo}
#' 
#' @examples \donttest{airq <- airquality[complete.cases(airquality), ]
#' set.seed(42)
#' airq.ens <- pre(Ozone ~ ., data = airquality[complete.cases(airquality),])
#' singleplot(airq.ens, "Temp")
#' 
#' ## For multinomial and mgaussian families, one PDP is created per category or outcome
#' set.seed(42)
#' airq.ens3 <- pre(Ozone + Wind ~ ., data = airq, family = "mgaussian")
#' singleplot(airq.ens3, varname = "Day")
#' 
#' set.seed(42)
#' iris.ens <- pre(Species ~ ., data = iris, family = "multinomial")
#' singleplot(iris.ens, varname = "Petal.Width")}
#' @seealso \code{\link{pre}}, \code{\link{pairplot}}
#' @export
singleplot <- function(object, varname, penalty.par.val = "lambda.1se",
                       nvals = NULL, type = "response", ylab = NULL, 
                       response = NULL,
                       gamma = NULL, newdata = NULL, xlab = NULL, rug = TRUE, ...)
{
  
  ## Check if newdata supplied matches original data
  if (!is.null(newdata)) {
    if (!all(object$x_names %in% colnames(newdata))) {
      stop("Newdata must contain all predictors used to fit original ensemble.")
    }
  }
  
  ## Check if proper object argument is specified
  if (!inherits(object, "pre")) {
    stop("Argument object should be an object of class 'pre'")
  }
  
  ## Check if proper varname argument is specified
  if (length(varname) != 1L || !is.character(varname)) {
    stop("Argument varname should be a character vector of length 1.")
  } else if (!(varname %in% object$x_names)) {
    varnames <- grep(varname, x = object$x_names, value = TRUE, fixed = TRUE)
    if (length(varnames > 0)) {
      stop("Argument varname should specify the variable name as specified in the model formula (e.g., ", paste0(paste0("'", varnames, "'"), collapse = " or "), ").")
    } else {
      stop("Argument varname should specify the name of a variable used to generate the ensemble.")
    }
  }
  
  ## Check if proper penalty.par.val argument is specified
  if (length(penalty.par.val) != 1L) {
    stop("Argument penalty.par.val should be a vector of length 1.")
  } else if (!( penalty.par.val %in% c("lambda.min", "lambda.1se") || 
                (is.numeric(penalty.par.val) && penalty.par.val >= 0))) {
    stop("Argument penalty.par.val should be equal to 'lambda.min', 'lambda.1se' or a numeric value >= 0")
  }
  
  ## check if proper gamma value is specified
  if (!is.null(gamma)) {
    if (!is.null(object$glmnet.fit$relaxed)) { 
      warning("A value for gamma was specified, but will be ignored because the rule ensemble was not fit using relax = TRUE.")
      gamma <- NULL
    }
  }
  if (!is.null(gamma)) {
    if (!(length(gamma) == 1L && gamma >= 0 && gamma <= 1)) {
      stop("Argument gamma has been supplied, but should be a single numeric value [0, 1].")
    }
    if (!gamma %in% object$glmnet.fit$relaxed$gamma) {
      stop("Specified gamma value should be one of: ", 
           paste(object$glmnet.fit$relaxed$gamma, sep = ", "),
           ".")
    }
  }
  
  ## Check if proper nvals argument is specified
  if (!is.null(nvals)) {
    if (length(nvals) != 1L || nvals != as.integer(nvals)) {
      stop("Argument nvals should be an integer vector of length 1.")
    } else if (is.factor(object$data[ , varname]) && !is.null(nvals)) {
      warning("Plot is requested for variable of class factor. Value specified for
              nvals will be ignored, all factor levels will be used.", immediate. = TRUE)
      nvals <- NULL
    }
  }
  
  ## Check if proper type argument is specified
  if (length(type) != 1L || !is.character(type)) {
    stop("Argument type should be a single character string.")
  }
  
  ## Generate expanded dataset
  if (is.null(newdata)) newdata <- object$data
  if (is.null(nvals)) {
    newx <- unique(newdata[ , varname])
  } else {
    newx <- seq(
      min(newdata[ , varname]), max(newdata[ , varname]), length = nvals)
  }
  exp_dataset <- if (is.null(newdata)) {
    newdata[rep(row.names(newdata), times = length(newx)), ]
  } else {
    newdata[rep(row.names(newdata), times = length(newx)), ]
  }
  exp_dataset[ , varname] <- rep(newx, each = nrow(newdata))
  
  ## get predictions
  if (is.null(gamma)) {
    exp_dataset$predy <- predict.pre(object, newdata = exp_dataset, type = type,
                                     penalty.par.val = penalty.par.val)
  } else {
    exp_dataset$predy <- predict.pre(object, newdata = exp_dataset, type = type,
                                     penalty.par.val = penalty.par.val, 
                                     gamma = gamma)
  }
  
  ## create plot
  if (object$family %in% c("multinomial", "mgaussian")) {
    resp_names <- if (is.null(response)) {
      colnames(exp_dataset$predy) 
    } else {
      colnames(exp_dataset$predy)[response]
    }
    pd <- list()
    for (resp_name in resp_names) {
      y_lab <- ifelse(is.null(ylab), resp_name, ylab)
      pd[[resp_name]] <- aggregate(exp_dataset$predy[ , resp_name], 
                                   by = exp_dataset[varname], 
                                   data = exp_dataset, FUN = mean)
      plot(pd[[resp_name]], type = "l", ylab = y_lab, 
           xlab = if (is.null(xlab)) varname else xlab, ...)
      if (rug) rug(pd[[resp_name]][[varname]])
    }
  } else {
    pd <- aggregate(exp_dataset$predy, by = exp_dataset[varname], 
                    data = exp_dataset, FUN = mean)
    y_lab <- ifelse(is.null(ylab), object$y_names, ylab) 
    plot(pd, type = "l", ylab = y_lab, 
         xlab = if (is.null(xlab)) varname else xlab, ...)
    if (rug) rug(pd[[varname]])
  }
  invisible(pd)
}





#' Create partial dependence plot for a pair of predictor variables in a prediction 
#' rule ensemble (pre)
#'
#' \code{pairplot} creates a partial dependence plot to assess the effects of a
#' pair of predictor variables on the predictions of the ensemble. Note that plotting 
#' partial dependence is computationally intensive. Computation time will increase 
#' fast with increasing numbers of observations and variables. For large 
#' datasets, package `plotmo` (Milborrow, 2019) provides more efficient functions 
#' for plotting partial dependence and also supports `pre` models. 
#'
#' @param object an object of class \code{\link{pre}}
#' @param gamma Mixing parameter for relaxed fits. See  
#' \code{\link[glmnet]{coef.cv.glmnet}}.
#' @param varnames character vector of length two. Currently, pairplots can only
#' be requested for non-nominal variables. If varnames specifies the name(s) of
#' variables of class \code{"factor"}, an error will be printed.
#' @inheritParams print.pre
#' @param type character string. Type of plot to be generated. 
#' \code{type = "heatmap"} yields a heatmap plot, \code{type = "contour"} yields 
#' a contour plot, \code{type = "both"} yields a heatmap plot with added contours,
#' \code{type = "perspective"} yields a three dimensional plot.
#' @param nvals optional numeric vector of length 2. For how many values of
#' x1 and x2 should partial dependence be plotted? If \code{NULL}, a grid of all possible
#' combinations of the observed values of the two predictor variables specified will be used 
#' (see details).
#' @param response numeric vector of length 1. Only relevant for multivariate gaussian 
#' and multinomial responses. If \code{NULL} (default), PDPs for all response 
#' variables or categories will be produced. A single integer can be specified, 
#' indicating for which response variable or category PDPs should be produced.   
#' @param pred.type character string. Type of prediction to be plotted on z-axis.
#' \code{pred.type = "response"} gives fitted values for continuous outputs and
#' fitted probabilities for nominal outputs. \code{pred.type = "link"} gives fitted
#' values for continuous outputs and linear predictor values for nominal outputs.
#' @param xlab character. Label to be printed on the x-axis. If \code{NULL},
#' the first elements of the supplied \code{varnames} will be printed on the x-axis.
#' @param ylab character. Label to be printed on the y-axis. If \code{NULL},
#' the second element of the supplied \code{varnames} will be printed on the y-axis.
#' @param newdata Optional \code{data.frame} in which to look for variables 
#' with which to predict. If \code{NULL}, the \code{data.frame} used to fit the 
#' original ensemble will be used.
#' @param main character vector. Title(s) for the plot. If \code{NULL}, the name 
#' of the response will be printed.
#' @param rug logical. Ignored if \code{type = "perspective"}. Should a rug be 
#' plotted on the x- and y-axes, representing the location of observed datapoint? 
#' Note that the rugs will only show where values have been observed, not their 
#' frequency/density.
#' @param ... Further arguments to be passed to \code{\link[graphics]{image}}, 
#' \code{\link[graphics]{contour}} or \code{\link[graphics]{persp}} (depending on
#' whether \code{type} is specified to be \code{"heatmap"}, \code{"contour"}, \code{"both"} 
#' or \code{"perspective"}).
#' 
#' @details 
#' Partial dependence functions are described in section 8.1 of Friedman & 
#' Popescu (2008).
#' 
#' By default, partial dependence will be plotted for each combination
#' of 20 values of the specified predictor variables. When \code{nvals = NULL} is
#' specified, a dependence plot will be created for every combination of the unique
#' observed values of the two specified predictor variables. If \code{NA} instead of
#' a numeric value is specified for one of the predictor variables, all observed
#' values for that variables will be used. Specifying \code{nvals = NULL} and 
#' \code{nvals = c(NA, NA)} will yield the exact same result.
#' 
#' High values, \code{NA} or \code{NULL} for \code{nvals} result in long 
#' computation times and possibly memory problems. Also, \code{\link{pre}} 
#' ensembles derived from training datasets that are very wide or long may 
#' result in long computation times and/or memory allocation errors. 
#' In such cases, reducing
#' the values supplied to \code{nvals} will reduce computation time and/or
#' memory allocation errors. 
#' 
#' When numeric value(s) are specified for \code{nvals}, values for the
#' minimum, maximum, and nvals - 2 intermediate values of the predictor variable
#' will be plotted. 
#' 
#' Alternatively, \code{newdata} can be specified to provide a different (smaller) 
#' set of observations to compute partial dependence over.
#' If \code{mi_pre} was used to derive the original rule ensemble, 
#' \code{newdata = "mean.mi"} can be specified. This 
#' will result in an average dataset being computed over the imputed datasets, 
#' which are then used to compute partial dependence functions. This greatly 
#' reduces the number of observations and thereby computation time.
#' 
#' If none of the variables specified with argument \code{varnames} was
#' selected for the final prediction rule ensemble, an error will be returned.
#' 
#' @examples \donttest{airq <- airquality[complete.cases(airquality),]
#' set.seed(42)
#' airq.ens <- pre(Ozone ~ ., data = airq)
#' pairplot(airq.ens, c("Temp", "Wind"))
#' 
#' ## For multinomial and mgaussian families, one PDP is created per category or outcome
#' set.seed(42)
#' airq.ens3 <- pre(Ozone + Wind ~ ., data = airq, family = "mgaussian")
#' pairplot(airq.ens3, varnames = c("Day", "Month"))
#' 
#' set.seed(42)
#' iris.ens <- pre(Species ~ ., data = iris, family = "multinomial")
#' pairplot(iris.ens, varname = c("Petal.Width", "Petal.Length"))}
#' @export
#' @references Friedman, J. H., & Popescu, B. E. (2008). Predictive learning 
#' via rule ensembles. \emph{The Annals of Applied Statistics, 2}(3), 916-954.
#' 
#' Milborrow, S. (2019). plotmo: Plot a model's residuals, response, and partial 
#' dependence plots. \url{https://CRAN.R-project.org/package=plotmo}
#' @import graphics
#' @seealso \code{\link{pre}}, \code{\link{singleplot}} 
#' @export
pairplot <- function(object, varnames, type = "both", gamma = NULL,
                     penalty.par.val = "lambda.1se", response = NULL,
                     nvals = c(20L, 20L), pred.type = "response", 
                     newdata = NULL, xlab = NULL, ylab = NULL, main = NULL, 
                     rug = TRUE, ...) {
  
  ## check if package interp is installed
  if (!(requireNamespace("interp"))) {
    stop("Function pairplot requires package interp, install package interp first.")
  }
  
  ## Check if proper object argument is specified
  if (!inherits(object, "pre")) {
    stop("Argument object should be an object of class 'pre'")
  }
  
  ## Check if newdata supplied matches original data
  if (!is.null(newdata)) {
    if (!all(object$x_names %in% colnames(newdata))) {
      stop("Newdata must contain all predictors used to fit original ensemble.")
    }
  }
  
  ## Check if proper varnames argument is specified
  if (length(varnames) != 2L || !is.character(varnames)) {
    stop("Argument varnames should be a character vector of length 2.")
  } else if (!(all(varnames %in% object$x_names))) {
    varname <- grep(varnames[1], x = object$x_names, value = TRUE, fixed = TRUE)
    varnames <- c(varname, grep(varnames[2], x = object$x_names, value = TRUE, fixed = TRUE))
    if (length(varnames > 0)) {
      stop("Argument varnames should specify the variable names as specified in the model formula (e.g., ", paste0(paste0("'", varnames, "'"), collapse = " and/or "), ").")
    } else {
      stop("Argument varnames should specify names of variables used to generate the ensemble.")
    }
  } else if (any(sapply(object$data[ , varnames], is.factor))) {
    stop("3D partial dependence plots are currently not supported for factors.")
  }
  
  ## Check if proper penalty.par.val argument is specified
  if (!(length(penalty.par.val) == 1)) {
    stop("Argument penalty.par.val should be a vector of length 1.")
  } else if (!(penalty.par.val %in% c("lambda.min", "lambda.1se") || 
               (is.numeric(penalty.par.val) && penalty.par.val >= 0))) {
    stop("Argument penalty.par.val should be equal to 'lambda.min', 'lambda.1se' or a numeric value >= 0")
  }
  
  ## check if proper gamma value is specified
  if (!is.null(gamma)) {
    if (!is.null(object$glmnet.fit$relaxed)) { 
      warning("A value for gamma was specified, but will be ignored because the rule ensemble was not fit using relax = TRUE.")
      gamma <- NULL
    }
  }
  if (!is.null(gamma)) {
    if (!(length(gamma) == 1L && gamma >= 0 && gamma <= 1)) {
      stop("Argument gamma has been supplied, but should be a single numeric value [0, 1].")
    }
    if (!gamma %in% object$glmnet.fit$relaxed$gamma) {
      stop("Specified gamma value should be one of: ", 
           paste(object$glmnet.fit$relaxed$gamma, sep = ", "),
           ".")
    }
  }
  
  ## Check if proper nvals argument is specified 
  if (is.null(nvals)) nvals <- c(NA, NA)
  
  if (!any(is.na(nvals))) {
    if (!(length(nvals) == 2 && all(nvals == as.integer(nvals)))) {
      stop("Argument nvals should be an integer vector of length 2.")
    }
    if (length(unique(object$data[ , varnames[1]])) < nvals[1] ||
        length(unique(object$data[ , varnames[2]])) < nvals[2]) {
      uniq1 <- unique(object$data[ , varnames[1]])
      uniq2 <- unique(object$data[ , varnames[2]])
      warning(paste0("The nvals argument specified more values than the number of observed unique values of ",
                     varnames[1], " and/or ", varnames[2], 
                     ". Specifying nvals=NULL, or NA for one of the predictors may reduce computation time. "))
    }
  }
  
  ## Check if proper type argument is specified
  if (!(length(type) == 1 && is.character(type))) {
    stop("Argument type should be equal to 'heatmap', 'contour', 'both' or 'perspective'.")
  }
  
  ## Check if proper pred.type argument is specified
  if (!(length(type) == 1 && is.character(type))) {
    stop("Argument type should be a single character string.")
  }
  
  ## generate expanded dataset
  if (is.null(newdata)) newdata <- object$data
  if (is.null(nvals)) {
    newx1 <- sort(unique(newdata[ , varnames[1]]))
    newx2 <- sort(unique(newdata[ , varnames[2]]))
  } else {
    newx1 <- if(is.na(nvals)[1]) {
      sort(unique(unique(newdata[ , varnames[1]])))
    } else {
      seq(min(newdata[ , varnames[1]]), max(newdata[ , varnames[1]]),
          length = nvals[1])
    }
    newx2 <- if(is.na(nvals)[2]) {
      sort(unique(unique(newdata[ , varnames[2]])))
    } else {
      seq(min(newdata[ , varnames[2]]), max(newdata[ , varnames[2]]),
          length = nvals[2])
    }
  }
  nobs1 <- length(newx1)
  nobs2 <- length(newx2)
  nobs <- nobs1*nobs2
  exp_dataset <- newdata[rep(row.names(newdata), times = nobs), ]
  exp_dataset[ , varnames[1]] <- rep(newx1, each = nrow(newdata)*nobs2)
  exp_dataset[ , varnames[2]] <- rep(rep(newx2, each = nrow(newdata)),
                                     times = nobs1)
  
  ## compute predictions
  if (is.null(gamma)) {
    pred_vals <- predict.pre(object, newdata = exp_dataset, type = pred.type,
                             penalty.par.val = penalty.par.val)
  } else {
    pred_vals <- predict.pre(object, newdata = exp_dataset, type = pred.type,
                             penalty.par.val = penalty.par.val, gamma = gamma)
  }
  
  ## create plot
  if (object$family %in% c("mgaussian", "multinomial")) {
    resp_names <- if (is.null(response)) {
      colnames(pred_vals) 
    } else {
      colnames(pred_vals)[response]
    }
    pd <- list()
    if (!is.null(main)) {
      if (length(main) == 1L) main <- rep(main, times = length(resp_names))
    }
    for (resp_name in resp_names) {
      main_tmp <- ifelse(is.null(main), resp_name, main[which(resp_names) == resp_name])
      xyz <- interp::interp(exp_dataset[ , varnames[1]], exp_dataset[ , varnames[2]],
                            pred_vals[, resp_name], duplicate = "mean")
      pd[[resp_name]] <- xyz
      if (type == "heatmap" || type == "both") {
        if (is.null(match.call()$col)) {
          colors <- rev(c("#D33F6A", "#D95260", "#DE6355", "#E27449", "#E6833D", 
                          "#E89331", "#E9A229", "#EAB12A", "#E9C037", "#E7CE4C", 
                          "#E4DC68", "#E2E6BD"))
          image(xyz, xlab = if (is.null(xlab)) varnames[1] else xlab, 
                ylab = if (is.null(ylab)) varnames[2] else ylab, 
                col = colors, main = main_tmp, ...)
          if (rug) {
            rug(object$data[[varnames[1]]], side = 1)
            rug(object$data[[varnames[2]]], side = 2)
          }
        } else {
          image(xyz, xlab = if (is.null(xlab)) varnames[1] else xlab, 
                ylab = if (is.null(ylab)) varnames[2] else ylab, main = main_tmp, ...)
          if (rug) {
            rug(object$data[[varnames[1]]], side = 1)
            rug(object$data[[varnames[2]]], side = 2)
          }
        }
        if (type == "both") {
          contour(xyz, add = TRUE)
        }
      }
      if (type == "contour") {
        contour(xyz, xlab = if (is.null(xlab)) varnames[1L] else xlab, 
                ylab = if (is.null(ylab)) varnames[2L] else ylab, main = main_tmp, ...) 
        if (rug) {
          rug(object$data[[varnames[1]]], side = 1L)
          rug(object$data[[varnames[2]]], side = 2L)
        }
      }
      if (type == "perspective") {
        persp(xyz, xlab = if (is.null(xlab)) varnames[1] else xlab, 
              ylab = if (is.null(ylab)) varnames[2] else ylab, zlab = "predicted", 
              main = main_tmp, ...)
      }
    }
    invisible(pd)
  } else { ## family is not multinomial or mgaussian
    xyz <- interp::interp(exp_dataset[ , varnames[1]], exp_dataset[ , varnames[2]],
                          pred_vals, duplicate = "mean")
    main <- ifelse(is.null(main), object$y_names, main)
    if (type == "heatmap" || type == "both") {
      if (is.null(match.call()$col)) {
        colors <- rev(c("#D33F6A", "#D95260", "#DE6355", "#E27449", "#E6833D", 
                        "#E89331", "#E9A229", "#EAB12A", "#E9C037", "#E7CE4C", 
                        "#E4DC68", "#E2E6BD"))
        image(xyz, xlab = if (is.null(xlab)) varnames[1] else xlab, 
              ylab = if (is.null(ylab)) varnames[2] else ylab, 
              col = colors, main = main, ...)
        if (rug) {
          rug(object$data[[varnames[1L]]], side = 1L)
          rug(object$data[[varnames[2L]]], side = 2L)
        }
      } else {
        image(xyz, xlab = if (is.null(xlab)) varnames[1L] else xlab, 
              ylab = if (is.null(ylab)) varnames[2L] else ylab, main = main, ...)
        if (rug) {
          rug(object$data[[varnames[1L]]], side = 1L)
          rug(object$data[[varnames[2L]]], side = 2L)
        }
      }
      if (type == "both") {
        contour(xyz, add = TRUE)
      }
    }
    if (type == "contour") {
      contour(xyz, xlab = if (is.null(xlab)) varnames[1] else xlab, 
              ylab = if (is.null(ylab)) varnames[2] else ylab, main = main, ...) 
      if (rug) {
        rug(object$data[[varnames[1L]]], side = 1L)
        rug(object$data[[varnames[2L]]], side = 2L)
      }
    }
    if (type == "perspective") {
      persp(xyz, xlab = if (is.null(xlab)) varnames[1] else xlab, 
            ylab = if (is.null(ylab)) varnames[2] else ylab, zlab = "predicted",
            main = main, ...)
    }
    invisible(xyz)
  }
}