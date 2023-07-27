#' Fit a prediction rule ensemble to multiply-imputed data (experimental)
#'
#' Function \code{mi_pre} derives a sparse ensemble of rules and/or
#' linear rules based on imputed data. The function is still experimental,
#' so use at own risk.
#'
#' @inheritParams pre
#' @param data A list of imputed datasets. The datasets must have identically-named
#' columns, but need not have the same number of rows (this can happen, for example.
#' if a bootstrap sampling approach had been employed for multiple imputation).
#' @param weights A list of observation weights for each observation in each 
#' imputed dataset. The list must have the same length as \code{data}, and each 
#' element must be a numeric vector of length identical to the number of rows of 
#' the corresponding imputed dataset in \code{data}. The default is 
#' \code{NULL}, yielding constant observation weights w_i = 1/M, where M is the 
#' number of imputed datasets (i.e., \code{length(data)}).
#' @param obs_ids A list of observation ids, corresponding to the id in the
#' original data, of each observation in each imputed dataset. Defaults to 
#' \code{NULL}, which assumes that the imputed datasets contain the observations 
#' in identical order. If specified, the list must have
#' the same length as \code{data}, and each element must be a numeric or character 
#' vector of length identical to the number of rows of the corresponding imputed
#' dataset in \code{data}. At least some of the observations ids must be repeated 
#' at least some times, within or between imputed datasets.
#' @param compl_frac An optional list specifying the fraction of observed values
#' for each observation. This will be used to compute observation weights as
#' a function of the fraction of complete data per observations, as per 
#' Wan et al. (2015), but note that this is only recommended for users who
#' know the risks (i.e., an analysis more like complete-case analysis).
#' If specified, the list must have
#' the same length as \code{data}, and each element must be a numeric  
#' vector of length identical to the number of rows of the corresponding imputed
#' dataset in \code{data}. 
#'
#' @details Experimental function to fit a prediction rule ensemble to 
#' multiply imputed data. Essentially, it is a wrapper function around function
#' \code{pre()}, the main differences relate to sampling for the tree induction
#' and fold assignment for estimation of the coefficients for the final ensemble.
#' 
#' Function \code{mi_pre} implements a so-called stacking approach to the analysis
#' of imputed data (see also Wood et al., 2008), where imputed datasets are combined 
#' into one large dataset.
#' In addition to adjustments of the sampling procedures, adjustments to observation 
#' weight are made to counter the artificial inflation of sample size.  
#' 
#' Observations which occur repeatedly across the imputed datasets will be 
#' completely in- or excluded from each sample or fold, to avoid overfitting. Thus,
#' complete observations instead of individual imputed observations are sampled,
#' for tree and rule induction, as well as the cross-validation for selecting the
#' penalty parameter values for the final ensemble.
#' 
#' It is assumed that data have already been imputed (using e.g.,
#' R package mice or missForest), and therefore function \code{mi_pre} takes a 
#' \code{list} of imputed datasets as input data.
#' 
#' Although the option to use the fraction of complete data for computing 
#' observation weight is provided through argument \code{compl_frac}, users
#' are not advised to use it. See e.g., Du et al. (2022): "An alternative weight 
#' specification, proposed in Wan et al. (2015), is o_i = f_i/D, where f_i is 
#' the number of observed predictors out of the total number of predictors for 
#' subject i [...] upweighting subjects with less missingness and downweighting 
#' subjects with more missingness can, in some sense, be viewed as making the 
#' optimization more like complete-case analysis, which might be problematic 
#' for Missing at Random (MAR) and Missing not at Random (MNAR) scenarios."
#' 
#' @return An object of class \code{pre}.
#' @export
#' @seealso \code{\link{pre}} \code{\link{mi_mean}}
#' 
#' @references
#' Du, J., Boss, J., Han, P., Beesley, L. J., Kleinsasser, M., Goutman, S.A., ... 
#' & Mukherjee, B. (2022). Variable selection with multiply-imputed datasets: 
#' choosing between stacked and grouped methods. Journal of Computational and 
#' Graphical Statistics, 31(4), 1063-1075. \doi{10.1080/10618600.2022.2035739}.
#' 
#' Wood, A. M., White, I. R., & Royston, P. (2008). How should variable selection 
#' be performed with multiply imputed data? Statistics in medicine, 27(17), 
#' 3227-3246. \doi{10.1002/sim.3177}
#' 
#' @examples \donttest{library("mice")
#' set.seed(42)
#' 
#' ## Shoot extra holes in airquality data
#' airq <- sapply(airquality, function(x) {
#'   x[sample(1:nrow(airquality), size = 25)] <- NA
#'   return(x)
#' })
#' 
#' ## impute the data
#' imp <- mice(airq, m = 5)
#' imp <- as.list(complete(imp, action = "all"))
#' 
#' ## fit a rule ensemble to the imputed data
#' set.seed(42)
#' airq.ens.mi <- mi_pre(Ozone ~ . , data = imp)}
mi_pre <- function(formula, data, ## As in pre()
                   weights = NULL, 
                   obs_ids = NULL,
                   compl_frac = NULL,
                   nfolds = 10L,
                   sampfrac = .5,
                   ...) {
  
  if (!is.list(data))
    stop("Argument data must specify a list of imputed datasets.")
  
  if (is.null(obs_ids)) obs_ids <- lapply(data, rownames)
  if (!is.list(obs_ids)) {
    stop("Argument obs_ids must be a list.")
  } else if (length(obs_ids) != length(data)) {
    stop("Argument obs_ids should specify a list with the same number of elements as data.")
  } else if (!all(sapply(obs_ids, length) == sapply(data, nrow))) {
    stop("Each element of obs_ids should have length equal to the number of rows of the corresponding element of data.")
  } else if (length(unique(unlist(obs_ids))) == sum(sapply(data, nrow))) {
    warning("None of the observation ids specified by argument obs_ids seem to occur repeatedly.")
  }
  
  if (!is.null(compl_frac)) {  
    if (!is.list(compl_frac)) {
      stop("Argument compl_frac must be a list.")
    } else if (length(compl_frac) != length(data)) {
      stop("Argument compl_frac should specify a list with the same number of elements as data.")
    } else if (!all(sapply(compl_frac, length) == sapply(data, nrow))) {
      stop("Each element of compl_frac should have length equal to the number of rows of the corresponding element of data.")
    }
  }
  
  if (!all(unlist(unique(lapply(data, names))) == names(data[[1L]])))
    stop("All imputed datasets must contain the same variables.")
  
  M <- length(data)
  
  ## Process data
  data <- do.call(rbind, data)
  obs_ids <- unlist(obs_ids)
  
  ## Compute weights
  weights <- if (is.null(weights)) {
    ## Use constant weights
    rep(1/M, times = nrow(data))
  } else if (!is.null(compl_frac)) {
    unlist(compl_frac) / M
  } else {
    rep_len(weights, length.out = nrow(data))
  }
  
  ## Set up sampling function for tree induction
  samp_func <- function(n, weights) {
    samp.size <- round(sampfrac * length(unique(obs_ids))) 
    sampled.ids <- sample(unique(obs_ids), size = samp.size)
    which(obs_ids %in% sampled.ids)
  }
  
  ## Set up foldids for cv.glmnet
  ##
  ## Each unique observation id must first be assigned to a fold;
  ## First, foldid is a vector of length length(unique(unlist(obs_ids)));
  ## Next, copy the foldids to a vector of length nrow(data).
  ##
  foldid <- unique(obs_ids)
  names(foldid) <- foldid
  foldid[names(foldid)] <- sample(rep(1:nfolds, length = length(foldid)))
  foldid <- as.integer(foldid[obs_ids])
  
  ## Apply pre
  pre(formula = formula, data = data, weights = weights, 
      foldid = foldid, sampfrac = samp_func, ...)
}



#' Compute the average dataset over imputed datasets.
#' 
#' \code{mi_mean} computes the averages dataset over a list of imputed datasets.
#' Can be used to reduce computation time of functions \code{singleplot} and 
#' \code{pairplot}.
#' 
#' @param data List of imputed datasets to compute the average dataset over.
#' 
#' @details It is assumed every imputed dataset contains the same observations 
#' (but not the same values) in the same order.
#' @return A dataset that is the average over the imputed datasets specified
#' with \code{data}.
#' @export
#' @seealso \code{\link{mi_pre}}, \code{\link{singleplot}}, \code{\link{pairplot}}
mi_mean <- function(data) {
  data.frame(apply(simplify2array(lapply(data, as.matrix)), 1:2, mean))
}