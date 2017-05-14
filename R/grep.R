#' @export
gpre_tress <- function(
  ...,
  remove_duplicates_complements = TRUE,
  mtry = Inf, ntrees = 500,
  maxdepth = 3L, learnrate = 0.01,
  parallel = FALSE){
  if(learnrate < 0 && learnrate > 1)
    stop("learnrate must be between 0 and 1")
  
  if(learnrate > 0 && parallel)
    warning("Parallel will not be used with learnrate > 0 in gpre_tress")
  
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
          dat[subsample, ], response, control, ytrafo))
        # Collect rules from tree:
        rules <- c(rules, list.rules(tree))
      }
    } else {
      rules <- c()
      if(family == "gaussian"){
        mf <- model.frame(update(formula, . ~ -1), data = data)
        y_learn <- model.response(mf)
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
            dat[subsample, ], response, control, ytrafo))
          # Collect rules from tree:
          rules <- c(rules, list.rules(tree))
          # Substract predictions from current y:
          y_learn <- y_learn - learnrate * predict_party_minimal(
            tree, newdata = data)
        }
      } else if (family == "binomial"){
        data2 <- data.frame(data, offset = 0)
        stop("TODO: implement")
        glmtreeformula <-formula(
          paste(paste(y_name, " ~ 1 |"),
                paste(x_names, collapse = "+")))
        
        for(i in 1:ntrees) {
          # Take subsample of dataset:
          subsample <- sample_func(n = n, weights = weights)
          subsampledata <- data2[subsample,]
          # Grow tree on subsample:
          tree <- glmtree(glmtreeformula, data = subsampledata, family = "binomial", 
                          maxdepth = maxdepth + 1,  
                          offset = offset)
          # Collect rules from tree:
          rules <- c(rules, list.rules(tree))
          # Update offset:
          data2$offset <- data2$offset + learnrate * predict(
            tree, newdata = data2, type = "link")
        }
      } else 
        stop("family '", family, "' is not implemented for gpre_tress")
    }
    
    ###################
    ## Rules cleanup ##
    ###################
  
    rules <- unique(rules[rules != ""])
    # TODO: how does this work with factors?
    rules <- sort(unname(rules))
    rules <- paste0("rTerm(", rules, ")")
    
    if(remove_duplicates_complements){
      frm <- paste("~", paste0(rules, collapse = " + "))
      rulevars <- model.frame(stats::formula(frm), data)
      rulevars <- as.matrix(rulevars)
      
      # Remove duplicates
      duplicates <- which(duplicated(rulevars, MARGIN = 2))
      rulevars <- rulevars[, -duplicates]
      rules <- rules[-duplicates]
      
      # Remove compliments
      sds <- apply(rulevars, 2, sd)
      sds_distinct <- 
        sapply(unique(sds), function(x) c(x, sum(sds == x)))
      
      complements <- vector(mode = "logical", length(sds))
      for(i in 1:ncol(sds_distinct)){
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
    
    # TODO: how does this work with factors?
    c("-1", rules) 
  }
  
  out
}

#' @export
rTerm <- function(x){
  if(!is.logical(x))
    stop("Non-logical input passed to rule")
  
  attr(x, "description") <- deparse(substitute(x))
  x <- as.integer(x)
  class(x) <- "rTerm"
  x
}

#' @export
gpre_linear <- function(
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
      stop("Terms with higher order is not implemented in with gpre_linear")
    
    is_numeric_term <- attr(mt, "dataClasses")== "numeric"
    if(attr(mt, "response") > 0)
      is_numeric_term <- is_numeric_term & !seq_along(is_numeric_term) %in% attr(mt, "response")
    is_numeric_term <- which(is_numeric_term)    
    
    ####################################
    ## Winsorize if needed and return ##
    ####################################
    
    if(winsfrac == 0){
      return(paste0("lTerm(", names(is_numeric_term), ")"))
    }
    
    out <- sapply(is_numeric_term, function(i) {
      x <- mf[, i]
      x_name <- colnames(mf)[i]
      qs <- quantile(x, c(winsfrac, 1 - winsfrac))
      
      if(!normalize)
        return(
          paste0("lTerm(", x_name, 
                 ", lb = ", signif(qs[1], 2), 
                 ", ub = ", signif(qs[2], 2)))
      
      
      sd <- sd(pmax(pmin(x, qs[2]), qs[1]))
      paste0("lTerm(", x_name, 
             ", lb = ", signif(qs[1], 2), 
             ", ub = ", signif(qs[2], 2), 
             ", scale = ", signif(sd, 2), ")")
    })
    
    out
  }
}

#' @export
lTerm <- function(x, lb = -Inf, ub = Inf, scale = 1){
  if(!is.numeric(x))
    stop("lTerm must numeric")
  
  attr(x, "description") <- deparse(substitute(x))
  
  # The (arbitrary?) 0.4 is from
  # PREDICTIVE LEARNING VIA RULE ENSEMBLES
  x <- pmin(pmax(x, lb), ub) / scale * 0.4
  class(x) <- "lTerm"
  x
}

#' @export
gpre_earth <- function(
  ..., degree = 3, nk = 11, standardize = TRUE, 
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
      
      ts <- row.names(fit$cuts)[-1]
      ts <- gsub("h\\(", "\\(", ts)
      
      ts[interaction_degree == 1] <- 
        gsub("(\\(.+)\\)", "pmax\\1, 0)", ts[interaction_degree == 1])
      
      ts[interaction_degree > 1] <- 
        gsub("(\\([^\\)]+)\\)", "pmax\\1, 0)", ts[interaction_degree > 1])
      
      if(standardize){
        vars <- with(data, sapply(ts, function(x) eval(parse(text = x))))
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

#' @export
eTerm <- function(x, scale = 1){
  if(!is.numeric(x))
    stop("eTerm must numeric")
  
  attr(x, "description") <- deparse(substitute(x))
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
  out <- c(defaults, args[not_match])
  out$x <- x
  out$y <- y
  out$family <- family
  
  out
}

#' @export
gpre_sample <- function(sampfrac = .5){
  if(sampfrac <= 0 || sampfrac > 1)
    stop("sampfrac should be greater > 0 and <= 1")
  
  if(sampfrac == 1){
    return(function(n, weights){
      sample(1:n, size = n, replace = TRUE, prob = weights)
    })
  } else {
    return(function(n, weights){
      sample(1:n, size = round(sampfrac * n), replace = FALSE, prob = weights)
    })
  }
}

#' @export
gpre <- function(
  formula, data, 
  #type = "both", 
  base_learners = list(gpre_tress(), gpre_linear()),
  weights = rep(1, times = nrow(data)), 
  sample_func = gpre_sample(), learnrate = NULL,
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
      stop("gpre is only implemented for 2 levels factors")
    
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
  
  ##################################################
  ## Perform penalized regression on the ensemble ##
  ##################################################
  
  call_args <- get_cv.glmnet_args(
    args = cv.glmnet_args, x = x, y = y, family = family)
  
  glmnet.fit <- do.call(cv.glmnet, call_args)
  
  ####################
  ## Return results ##
  ####################
  
  result <- list(
    glmnet.fit = glmnet.fit, call = match.call, 
    family = family, base_learners = base_learners, 
    modmat_formula = glmnet_formula)
  
  if(model){
    result <- c(result, list(
      data = data, weights = weights))
  } else {
    result <- c(result, list(
      data = NULL, weights = NULL)) 
  }
  
  class(result) <- "gpre"
  result
}

#' @rdname print.pre
#' @export
#' @method print gpre
print.gpre <- function(x, penalty.par.val = "lambda.1se", ...){
  print.pre(x, penalty.par.val, ...)
}

#' @rdname coef.pre
#' @export
#' @method coef gpre
coef.gpre <- function(object, penalty.par.val = "lambda.1se", ...)
{
  coefs <- as(coef.glmnet(object$glmnet.fit, s = penalty.par.val, ...), 
              Class = "matrix")
  colnames(coefs) <- "coefficient"
  coefs <- coefs[coefs != 0, ,drop = FALSE]
  data.frame(description = row.names(coefs), coefs)
}