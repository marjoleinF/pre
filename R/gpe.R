#' @title Learner Functions Generators for gpe
#' 
#' @description 
#' Functions to get "learners" functions for \code{\link{gpe}}.
#' 
#' @param ... Unused
#' @param remove_duplicates_complements \code{TRUE} if rules with complementary or duplicate support should be removed from trees
#' @param mtry Depth of trees. The argument is passed the tree methods in \code{partykit} package
#' @param ntrees Number of trees to fit
#' @param maxdepth Maximum dpeth of trees 
#' @param learnrate Learning rate for methods. It is the \eqn{\nu} parameter in Friedman & Popescu (2008)
#' @param parallel \code{TRUE} if basis functions should be found in parallel
#' @param use_L2_loss \code{TRUE} if binary outcomes should use L2 loss instead logistic loss when \code{learnrate > 0}. That is, use \code{\link{ctree}} instead of \code{\link{glmtree}}
#' @param winsfrac Quantiles to winsorize linear terms at. The value should be in \eqn{[0,0.5)}
#' @param normalize \code{TRUE} if value should be scaled by \eqn{0.4} times the inverse standard deviation. The \eqn{0.4} is from Friedman & Popescu (2008)
#' @param degree Maximum degree of interactions in \code{\link{earth}} model
#' @param nk Maximum number of basis functions in \code{\link{earth}} model
#' @param ntrain Number of models to fit 
#' 
#' @details 
#' \code{gpe_tress} provides learners for tree method. Either \code{\link{ctree}} or \code{\link{glmtree}} from the \code{partykit} package will be used.
#' 
#' \code{gpe_linear} provides linear terms for the \code{gpe}.
#' 
#' \code{gpe_earth} provides basis functions where each factor is a hinge function. The model is estimated with \code{\link{earth}}.
#' 
#' @return 
#' A function that has formal arguments \code{formula, data, weights, sample_func, verbose, family, ...}. The function returns a vector with character where each element is a term for the final formula in the call to \code{\link{cv.glmnet}}
#' 
#' @seealso 
#' \code{\link{gpe}}, \code{\link{rTerm}}, \code{\link{lTerm}}, \code{\link{eTerm}}
#' 
#' @references 
#' Hothorn, T., & Zeileis, A. (2015). partykit: A modular toolkit for recursive partytioning in R. \emph{Journal of Machine Learning Research}, 16, 3905-3909.
#' 
#' Friedman, J. H. (1991). Multivariate adaptive regression splines. \emph{The annals of statistics}, 1-67.
#' 
#' Stanford University. Laboratory for Computational Statistics, & Friedman, J. H. (1993). Fast MARS.
#' 
#' Friedman, J. H., & Popescu, B. E. (2008). Predictive learning via rule ensembles. \emph{The Annals of Applied Statistics}, 916-954.
#' 
#' @export
gpe_tress <- function(
  ...,
  remove_duplicates_complements = TRUE,
  mtry = Inf, ntrees = 500,
  maxdepth = 3L, learnrate = 0.01,
  parallel = FALSE, use_L2_loss = FALSE){
  if(learnrate < 0 && learnrate > 1)
    stop("learnrate must be between 0 and 1")
  
  if(learnrate > 0 && parallel)
    warning("Parallel will not be used with learnrate > 0 in gpe_tress")
  
  out <- function(
    formula, data, weights, sample_func, verbose, family, ...){
    ################
    ## Find rules ##
    ################
  
    if(learnrate == 0) { # always use ctree()
      if(parallel)
        stop("Not implemented")
      
      input <- ctree_setup(formula, data = data, maxdepth = maxdepth, mtry = mtry)
      rules <- c()
      n <- nrow(data)
      for(i in 1:ntrees) {
        # Take subsample of dataset
        subsample <- sample_func(n = n, weights = weights)
        # Grow tree on subsample:
        #tree <- ctree(formula, data = data[subsample,], maxdepth = maxdepth, 
        #                mtry = mtry)
        tree <- with(input, ctree_minmal(
          dat[subsample, ], response, control, ytrafo, terms))
        # Collect rules from tree:
        rules <- c(rules, list.rules(tree))
      }
    } else {
      rules <- c()
      if(family == "gaussian" || (
        family == "binomial" && use_L2_loss)){
        mf <- model.frame(update(formula, . ~ -1), data = data)
        y_learn <- model.response(mf)
        if(family == "binomial"){
          if(length(levels(y_learn)) != 2)
            stop("Factor for outcome in must have two levels in gpe_tress with a learning rate")
          
          y_learn <- as.numeric(y_learn == levels(y_learn)[1])
          mt <- terms(mf)
          data[, as.character(attr(mt,"variables")[[2]])] <- y_learn
        }
        input <- ctree_setup(formula, data = data, maxdepth = maxdepth, mtry = mtry)
        n <- nrow(data)
        
        for(i in 1:ntrees) {
          # Take subsample of dataset
          subsample <- sample_func(n = n, weights = weights)
          # Grow tree on subsample:
          #tree <- ctree(formula, data = data[subsample,], maxdepth = maxdepth, 
          #                mtry = mtry)
          input$dat[subsample, input$response] <- y_learn[subsample]
          tree <- with(input, ctree_minmal(
            dat[subsample, ], response, control, ytrafo, terms))
          # Collect rules from tree:
          rules <- c(rules, list.rules(tree))
          # Substract predictions from current y:
          y_learn <- y_learn - learnrate * predict_party_minimal(
            tree, newdata = data)
        }
      } else if (family == "binomial"){
        data2 <- data.frame(data, offset = 0)
        mt <- terms(formula, data = data)
        
        if(attr(mt, "response") != 1)
          stop("Left hand site of formula must have one term")
        
        glmtreeformula <-stats::formula(paste0(
          as.character((attr(mt, "variables")[[2]])), " ~ 1 |", paste0(
            attr(mt, "term.labels"), collapse = " + ")))
        n <- nrow(data)
        
        for(i in 1:ntrees) {
          # Take subsample of dataset:
          subsample <- sample_func(n = n, weights = weights)
          subsampledata <- data2[subsample,]
          # Grow tree on subsample:
          tree <- glmtree(
            glmtreeformula, data = subsampledata, family = "binomial", 
            maxdepth = maxdepth + 1,  
            offset = offset, 
            epsilon = 1e-4) # we set the relative change in the deviance lower
                            # to speed up the computations
          # Collect rules from tree:
          rules <- c(rules, list.rules(tree))
          # Update offset:
          data2$offset <- data2$offset + learnrate * predict(
            tree, newdata = data2, type = "link")
        }
      } else 
        stop("family '", family, "' is not implemented for gpe_tress")
    }
    
    ###################
    ## Rules cleanup ##
    ###################
  
    rules <- unique(rules[rules != ""])
    rules <- sort(unname(rules))
    rules <- paste0("rTerm(", rules, ")")
    
    if(remove_duplicates_complements){
      frm <- paste("~", paste0(rules, collapse = " + "))
      rulevars <- model.frame(stats::formula(frm), data)
      rulevars <- as.matrix(rulevars)
      
      # Remove duplicates
      duplicates <- which(duplicated(rulevars, MARGIN = 2))
      if(length(duplicates) > 0){
        rulevars <- rulevars[, -duplicates]
        rules <- rules[-duplicates]
      }
      
      # Remove compliments
      sds <- apply(rulevars, 2, sd)
      sds_distinct <- 
        sapply(unique(sds), function(x) c(x, sum(sds == x)))
      
      complements <- vector(mode = "logical", length(sds))
      for(i in seq_len(ncol(sds_distinct))){
        if(sds_distinct[2, i] < 2)
          next
        
        indices <- which(sds == sds_distinct[1, i])
        for(j in 2:length(indices)){
          indices_prev <- indices[1:(j - 1)] 
          complements[indices_prev] <- 
            complements[indices_prev] | apply(
              rulevars[, indices_prev, drop = F] != rulevars[, indices[j]], 2, all)
        }
      }
      
      rules <- rules[!complements]
    }
    
    c(rules) 
  }
  
  out
}

#' @title Wrapper Functions for Terms in gpe
#' 
#' @description 
#' Wrapper functions for terms in gpe.
#' 
#' @param x Input symbol 
#' @param lb Lower quantile when winsorizing. \code{-Inf} yields no winsorizing in the lower tail
#' @param ub Lower quantile when winsorizing. \code{Inf} yields no winsorizing in the upper tail
#' @param scale Inverse value to time \code{x} by. Usually the standard deviation is used. \eqn{0.4 / scale} is used as the multiplier as suggested in Friedman & Popescu (2008)
#' 
#' @details 
#' The motivation to use wrappers is to ease getting the different terms as shown in the examples and simplify the formula passed to \code{\link{cv.glmnet}} in \code{\link{gpe}}. \code{lTerm} potentially rescale and/or winsorize \code{x} depending on the input. \code{eTerm} potentially rescale \code{x} depending on the input.
#' 
#' @return 
#' \code{x} potentially transformed with additional information provided in the attributes.
#' 
#' @examples
#' mt <- terms(
#' ~ rTerm(x1 < 0) + rTerm(x2 > 0) + lTerm(x3) + eTerm(x4), 
#' specials = c("rTerm", "lTerm", "eTerm"))
#' attr(mt, "specials")
#' # $rTerm
#' # [1] 1 2
#' # 
#' # $lTerm
#' # [1] 3
#' # 
#' # $eTerm
#' # [1] 4
#' 
#' @references
#' 
#' Friedman, J. H., & Popescu, B. E. (2008). Predictive learning via rule ensembles. \emph{The Annals of Applied Statistics}, 916-954.
#' 
#' @seealso 
#' \code{\link{gpe}}, \code{\link{gpe_tress}} \code{\link{gpe_linear}} \code{\link{gpe_earth}}
#' 
#' @export
rTerm <- function(x){
  if(!is.logical(x))
    stop("Non-logical input passed to rule")
  
  attr(x, "description") <- deparse(substitute(x))
  x <- as.integer(x)
  class(x) <- "rTerm"
  x
}

#' @rdname gpe_tress
#' @export
gpe_linear <- function(
  ..., winsfrac = .025, normalize = TRUE){
  if(winsfrac < 0 && winsfrac > 1)
    stop("winsfrac must be 0 <= winsfrac <= 1")
  
  function(formula, data, weights, sample_func, verbose, family, ...){
    ########################
    ## Find numeric terms ##
    ########################
    
    mf <- model.frame(formula, data)
    mt <- attr(mf, "terms")
    
    if(any(attr(mf, "order") > 1))
      stop("Terms with higher order is not implemented in with gpe_linear")
    
    is_numeric_term <- attr(mt, "dataClasses")== "numeric"
    if(attr(mt, "response") > 0)
      is_numeric_term <- is_numeric_term & !seq_along(is_numeric_term) %in% attr(mt, "response")
    is_numeric_term <- which(is_numeric_term)    
    
    ####################################
    ## Winsorize if needed and return ##
    ####################################
    
    if(winsfrac == 0){
      if(!normalize)
        return(paste0("lTerm(", names(is_numeric_term), ")"))
      
      sds <- apply(mf[, names(is_numeric_term)], 2, sd)
      out <- mapply(function(x, s) paste0("lTerm(", x, ", scale = ", s, ")"), 
                    x = names(is_numeric_term), s = signif(sds, 2))
      return(out)
    }
    
    out <- sapply(is_numeric_term, function(i) {
      x <- mf[, i]
      x_name <- colnames(mf)[i]
      qs <- quantile(x, c(winsfrac, 1 - winsfrac))
      
      if(!normalize)
        return(
          paste0("lTerm(", x_name, 
                 ", lb = ", signif(qs[1], 2), 
                 ", ub = ", signif(qs[2], 2), ")"))
      
      
      sd <- sd(pmax(pmin(x, qs[2]), qs[1]))
      paste0("lTerm(", x_name, 
             ", lb = ", signif(qs[1], 2), 
             ", ub = ", signif(qs[2], 2), 
             ", scale = ", signif(sd, 2), ")")
    })
    
    out
  }
}

#' @rdname rTerm
#' @export
lTerm <- function(x, lb = -Inf, ub = Inf, scale = 1 / 0.4){
  if(!is.numeric(x))
    stop("lTerm must numeric")
  
  attr(x, "description") <- deparse(substitute(x))
  attr(x, "lb") <- lb
  attr(x, "ub") <- ub
  attr(x, "scale") <- scale
  
  # The (arbitrary?) 0.4 is from
  # PREDICTIVE LEARNING VIA RULE ENSEMBLES
  x <- pmin(pmax(x, lb), ub) / scale * 0.4
  class(x) <- "lTerm"
  x
}

#' @rdname gpe_tress
#' @export
gpe_earth <- function(
  ..., degree = 3, nk = 11, normalize = TRUE, 
  ntrain = 100, learnrate = 0.01){
  
  if(learnrate < 0 && learnrate > 1)
    stop("learnrate must be between 0 and 1")
  
  out <- function(formula, data, weights, sample_func, verbose, family, ...){
    ###########
    ## Setup ##
    ###########
    
    n <- nrow(data)
    mf <- model.frame(formula, data = data)
    mt <- attr(mf, "terms")
    x <- model.matrix(mt, mf)
    y <- model.response(mf)
    
    if(family == "binomial"){
      message("Beware that gpe_earth will use L2 loss to train")
      y <- as.numeric(y)
    }
    
    # We later need to take care of the factor terms
    factor_terms <- which(
      attr(mt, "dataClasses")[
        (1 + attr(mt, "response")):length(attr(mt, "dataClasses"))] == 
        "factor")
    n_factors <- length(factor_terms)
    factor_terms <- names(factor_terms)
    
    if(n_factors > 0){
      add_escapes <- function(regexp)
        stringr::str_replace_all(regexp, "(\\W)", "\\\\\\1")
      
      factor_terms_regexp <- add_escapes(factor_terms)
      factor_labels <- lapply(mf[, factor_terms, drop = FALSE], levels)
      
      regexp_find <- list()
      regexp_replace <- list()
      for(i in 1:n_factors){
        regexp_find[[i]] <- add_escapes(paste0(
          factor_terms[i], factor_labels[[i]]))
        regexp_replace[[i]] <- add_escapes(paste0(
          "(", factor_terms[i], " == '", factor_labels[[i]], "')"))
      }
    }
    
    basis_funcs <- c()
    for(i in 1:ntrain){
      ##########################
      ## Find basis functions ##
      ##########################
      
      subsample <- sample_func(n = n, weights = weights)
      
      fit <- earth(
        x = x[subsample, , drop = FALSE], y = y[subsample], degree = degree, 
        nk = nk, pmethod = "none")
      
      if(learnrate > 0)
        y <- drop(y - learnrate * predict(fit, type = "response", newdata = x))
      
      ###########################################
      ## Format basis functions terms & return ##
      ###########################################
      
      # For details on the earth object see ?earth.object. The two key elements
      # are dirs and cuts
      
      # -1 for the intercept
      interaction_degree <- rowSums(fit$dirs[-1, ] != 0)
      
      # Replace 
      #   h(xyz)
      # with 
      #   pmax(xyz, 0)
      ts <- row.names(fit$cuts)[-1]
      ts <- gsub("h(\\([^\\)]+)\\)($|\\*)", "pmax\\1, 0\\)\\2", ts)
      
      # Check if we have factor terms and adjust these. That is, we replace 
      #   facxyz
      # with 
      #   (fac == 'xyz')
      if(n_factors > 0){
        has_factor <- sapply(factor_terms_regexp, grepl, x = ts, perl = TRUE)
        if(is.vector(has_factor)) has_factor <- t(has_factor)
        
        if(any(has_factor)){
          for(i in 1:n_factors){
            needs_replace <- which(has_factor[, i])
            if(length(needs_replace) == 0)
              next
            
            r_find <- regexp_find[[i]]
            r_replace <- regexp_replace[[i]]
            
            for(j in seq_along(r_find))
              ts[needs_replace] <- stringr::str_replace(
                ts[needs_replace], r_find[j], r_replace[j])
          }
        }
      }
      
      if(normalize){
        vars <- with(data, eval(parse(
          text = paste0("cbind(", paste0(ts, collapse = ", "), ")"))))
        sds <- apply(vars, 2, sd)
        
        ts <- mapply(
          function(x, s) paste0("eTerm(", x, ", scale = ", s, ")"),
          x = ts, s = signif(sds, 2))
      } else {
        ts <- paste0("eTerm(", ts, ")")
      }
      
      basis_funcs <- c(basis_funcs, ts)
    }
    
    basis_funcs <- unique(basis_funcs)
    basis_funcs
  }
  
  out
}

#' @rdname rTerm
#' @export
eTerm <- function(x, scale = 1 / 0.4){
  if(!is.numeric(x) && !is.logical(x))
    stop("eTerm must numeric")
  
  attr(x, "description") <- deparse(substitute(x))
  attr(x, "scale") <- scale
  
  # The (arbitrary?) 0.4 is from
  # PREDICTIVE LEARNING VIA RULE ENSEMBLES
  x <- x / scale * 0.4
  class(x) <- "eTerm"
  x
}

get_cv.glmnet_args <- function(args, x, y, weights, family){
  defaults <- list(
    nfolds =  10L, standardize = FALSE, 
    type.measure = "deviance", thres = 1e-07, 
    parallel = FALSE)
  
  not_match <- !(names(args) %in% names(defaults))
  # TODO: matching elements should replace default
  out <- c(defaults, args[not_match])
  out$x <- x
  out$y <- y
  out$weights <- weights
  out$family <- family
  
  out
}

#' @title Sampling Function Generator for gpe
#' 
#' @description 
#' Provides a sample function for \code{\link{gpe}}.
#'
#' @param sampfrac Fraction of \code{n} to use for sampling. It is the \eqn{\eta / N} in Friedman & Popescu (2008)
#' 
#' @return 
#' Returns a function that takes an \code{n} argument for the number of observation and an \code{weights} argument for the case weights. The function return a vector of indices
#' 
#' @references
#' 
#' Friedman, J. H., & Popescu, B. E. (2008). Predictive learning via rule ensembles. \emph{The Annals of Applied Statistics}, 916-954.
#'
#' @seealso 
#' \code{\link{gpe}}
#' 
#' @export
gpe_sample <- function(sampfrac = .5){
  if(sampfrac <= 0 || sampfrac > 1)
    stop("sampfrac should be greater > 0 and <= 1")
  
  if(sampfrac == 1){
    return(function(n, weights){
      sample(1:n, size = n, replace = TRUE, prob = weights)
    })
  } else {
    has_written_about_weights <- FALSE
    return(function(n, weights){
      # Sub sampling will be used if all weights match
      all_weights_match <- all(weights[1] == weights)
      
      if(!all_weights_match && !has_written_about_weights){
        has_written_about_weights <<- TRUE
        message("Some weights do not match. Bootstrap will be used instead of subsampling to reflect weights")
      }
      
      sample(1:n, size = round(sampfrac * n), 
             replace = !all_weights_match, 
             prob = weights)
    })
  }
}

#' @title Derive a General Prediction Ensemble (gpe)
#' 
#' @description 
#' Provides an interface sparse prediction ensemble where basis functions are removed with the L1 penalty.
#' 
#' @param formula Symbolic description of the model to be fit of the form \code{y ~ x1 + x2 + ...+ xn}. If the output variable (left-hand side of the formula) is a factor, an ensemble for binary classification is created. Otherwise, an ensemble for prediction of a continuous variable is created
#' @param data \code{data.frame} containing the variables in the model
#' @param base_learners List of functions which has formal arguments \code{formula, data, weights, sample_func, verbose} and \code{family} and returns a vector of characters with terms for the final formula passed to \code{cv.glmnet}. See \code{\link{gpe_linear}}, \code{\link{gpe_tress}}, and \code{\link{gpe_earth}}
#' @param weights Case weights with length equal to number of rows in \code{data}
#' @param sample_func Function used to sample when learning with base learners. The function should have formal argument \code{n} and \code{weights} and return a vector of indices. See \code{\link{gpe_sample}}
#' @param verbose \code{TRUE} if comments should be posted throughout the computations
#' @param cv.glmnet_args List of arguments to \code{\link{cv.glmnet}}. \code{x}, \code{y}, \code{weights} and \code{family} will not be used
#' @param model \code{TRUE} if the \code{data} should added to the returned object
#' 
#' @details 
#' Provides are more general framework for making a sparse prediction ensemble than \code{\link{pre}}. A similar fit to \code{\link{pre}} can be estimated with the following call:
#' 
#' \code{
#' gpe(formula = y ~ x1 + x2 + x3, data = data, base_learners = list(gpe_linear(), gpe_tress()))
#'}
#'     
#' Products of hinge functions using MARS can be added to the ensemble above with the following calling:
#' 
#' \code{
#' gpe(formula = y ~ x1 + x2 + x3, data = data, base_learners = list(gpe_linear(), gpe_tress(), gpe_earth))
#' }
#' 
#' Other customs base learners can be implemented. See \code{\link{gpe_tress}} \code{\link{gpe_linear}} or \code{\link{gpe_earth}} for details of the setup. The sampling function given by \code{sample_func} can also be replaced by a custom sampling function. See \code{\link{gpe_sample}} for details of the setup.
#' 
#' @return 
#' An object of class \code{gpe}
#' 
#' @seealso 
#' \code{\link{pre}}, \code{\link{gpe_tress}}, \code{\link{gpe_linear}}, \code{\link{gpe_earth}}, \code{\link{gpe_sample}}
#' 
#' @references 
#' Friedman, J. H., & Popescu, B. E. (2008). Predictive learning via rule ensembles. \emph{The Annals of Applied Statistics}, 916-954.
#' 
#' @export
gpe <- function(
  formula, data, 
  base_learners = list(gpe_tress(), gpe_linear()),
  weights = rep(1, times = nrow(data)), 
  sample_func = gpe_sample(),
  verbose = FALSE, cv.glmnet_args = list(), 
  model = TRUE){
  
  ###################
  ## Preliminaries ##
  ###################
  
  if (!is.data.frame(data)) {
    stop("data should be a data frame.")
  }
  if (!is.logical(verbose)) {
    stop("Bad value for 'verbose'.")
  }
  
  mf <- model.frame(update(.~1, formula), data)
  y <- model.response(mf)
  n <- nrow(data)
  
  if (is.factor(y)) {
    if(length(levels(y)) != 2)
      stop("gpe is only implemented for 2 levels factors")
    
    family <- "binomial"
  } else {
    family <- "gaussian"
  }
  
  if(!all(unlist(lapply(base_learners, is.function))))
    stop("All the elements in base_learners must be functions")
  
  if(!is.function(sample_func))
    stop("sample_func must be a function")
  
  ############################
  ## Derive basis functions ##
  ############################
  
  formulas <- lapply(
    base_learners, function(f) 
      f(formula = formula, data = data, weights = weights,
        sample_func = sample_func, verbose = verbose, family = family))
  
  glmnet_formula <- lapply(formulas, paste0, collapse = " + ")
  glmnet_formula <- paste0(unlist(glmnet_formula), collapse = " + ")
  glmnet_formula <- stats::formula(paste("~", glmnet_formula))
  x <- model.matrix(glmnet_formula, data = data)
  Terms <- terms(glmnet_formula, data = data)
  
  ##################################################
  ## Perform penalized regression on the ensemble ##
  ##################################################
  
  call_args <- get_cv.glmnet_args(
    args = cv.glmnet_args, x = x, y = y, family = family, weights = weights)
  
  glmnet.fit <- do.call(cv.glmnet, call_args)
  
  ####################
  ## Return results ##
  ####################
  
  result <- list(
    glmnet.fit = glmnet.fit, call = match.call, 
    family = family, base_learners = base_learners, 
    modmat_formula = glmnet_formula, terms = Terms)
  
  if(model){
    result <- c(result, list(
      data = data, weights = weights))
  } else {
    result <- c(result, list(
      data = NULL, weights = NULL)) 
  }
  
  class(result) <- "gpe"
  result
}