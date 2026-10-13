##########
##
## Internal function to remove duplicate and complement rules:
##
delete_duplicates_complements <- function(
  rules, data, removecomplements = TRUE, removeduplicates = TRUE,
  return.dupl.compl = FALSE, sparse = FALSE, keep_rulevars = FALSE) {
  
  ## Generate rule variables:
  rulevars <- if(sparse)
    .get_rules_mat_sparse(data, rules) else
      .get_rules_mat_dense(data, rules)
  rules_to_process <- length(rules) > 0
  if (rules_to_process) {
    colnames(rulevars) <- names(rules) <- paste0("rule", 1:length(rules))
  }
  
  ## Remove duplicate rules
  if (removeduplicates && rules_to_process) {
    # Remove rules with identical support
    duplicates <- duplicated(rulevars, MARGIN = 2)
    duplicates.removed <- rules[duplicates]
    rulevars <- rulevars[, !duplicates, drop = FALSE]
    rules <- rules[!duplicates]
    
  } else 
    duplicates.removed <- NULL
  
  ## Remove complement rules
  if (removecomplements && rules_to_process) {
    if (sparse) {
      is_all_binary <- all(rulevars@x %in% 0:1)
      if (!is_all_binary) stop("method not implemented for non-binary rules")
      # find columns which length of non-zero entries are equal to the number 
      # of rows. The latter are complements if the union of the row indices are
      # equal to the number of rows
      js <- rep(1:ncol(rulevars), diff(rulevars@p))
      is <- rulevars@i + 1
      n <- rulevars@Dim[1]
      p <- rulevars@Dim[2]
      
      Jfac <- factor(js, levels = 1:p)
      row_indices <- split(is, Jfac)
      lengths <- table(Jfac)
      complements <- logical(p)
      for(i in 1:(p - 1)) {
        if (complements[i])
          next
        
        is_potential <- which(lengths[i] + lengths[(i + 1):p] == n) + i
        is_potential <- is_potential[!complements[is_potential]]
        
        if (length(is_potential) == 0)
          next
        
        union_len <- sapply(
          mapply(union, x = row_indices[i], 
                 y = row_indices[is_potential], SIMPLIFY = FALSE), 
          length)
        is_compl <- which(union_len == n) 
        if (length(is_compl) > 0) {
          complements[is_potential[is_compl]] <- TRUE
        }
      }
        
    } else { # sparse = FALSE
      if (!is.logical(rulevars)) {
        stop("method not implemented for non-binary rules")
      }
      
      # find columns will equal variance to reduce the number of comparisons
      vars <- apply(rulevars, 2, var_bin)
      vars_distinct <- lapply(
        unique(vars), function(x) { 
          idx <- which(is_almost_eq(x, vars))
          list(var = x, n = length(idx), idx = idx)
        })
      
      complements <- logical(ncol(rulevars))
      for (va in vars_distinct) {
        if(va$n < 2L)
          next
        
        idx <- va$idx
        idx <- setdiff(idx, which(complements))

        if(length(idx) < 2)
          next
        
        n_idx <- length(idx)
        for(j in 1:(n_idx - 1)){
          if (complements[idx[j]])
            next
          
          this_val <- rulevars[, idx[j]]
          is_compl <- 
            which(apply(
              rulevars[, idx[(j + 1):n_idx], drop = FALSE], 2, 
              function(x) all(x != this_val))) + j
          
          if (length(is_compl) > 0)
            complements[idx[is_compl]] <- TRUE
        }
      }
      
    }
    
    complements <- which(complements)
    complements.removed <- rules[complements]
    if (length(complements) > 0)
      rules <- rules[-complements]
    
  } else {
    complements.removed <- NULL
  }
  
  rulevars = if (keep_rulevars && rules_to_process)
    rulevars[, names(rules), drop = FALSE] else NULL
  
  ## Return results:
  if (return.dupl.compl) {
    return(list(
      rules = rules, rulevars = rulevars,
      duplicates.removed = duplicates.removed,
      complements.removed = complements.removed))
  }
  
  list(rules = rules, rulevars = rulevars)
}

# see https://stackoverflow.com/a/51457395/5861244
duplicated.dgCMatrix <- function (x, incomparables = NULL, MARGIN, ...) {
  ## incomparables = NULL added to match duplicated generic: function(x, incomparables, ...)
  MARGIN <- as.integer(MARGIN)
  n <- nrow(x)
  p <- ncol(x)
  J <- rep(1:p, diff(x@p))
  I <- x@i + 1
  X <- x@x
  if (MARGIN == 1L) {
    ## check duplicated rows
    names(X) <- J
    RowLst <- split(X, I)
    is_empty <- setdiff(1:n, I)
    result <- duplicated.default(RowLst)
  } else if (MARGIN == 2L) {
    ## check duplicated columns
    names(X) <- I
    ColLst <- split(X, J)
    is_empty <- setdiff(1:p, J)
    result <- duplicated.default(ColLst)
  } else {
    warning("invalid MARGIN; return NULL")
    result <- NULL
  }
  
  if(any(is_empty)){
    out <- logical(if(MARGIN == 1L) n else p)
    out[-is_empty] <- result
    if(length(is_empty) > 1)
      out[is_empty[-1]] <- TRUE
    result <- out
  }
  
  result
}

# see https://stackoverflow.com/a/30089750/5861244
cbind_sparse_vec <- function (...) {
  args <- list(...)
  stopifnot(all(sapply(args, inherits, what = "dsparseVector")))
  
  un_lengths <- unique(sapply(args,length))
  stopifnot(length(un_lengths) == 1)
  
  return(sparseMatrix( 
    x = unlist(lapply(args, slot, "x")), 
    i = unlist(lapply(args, slot ,"i")), 
    p = c(0, cumsum(sapply(args, function(x) { length(x@x) } ))),
    dims=c(un_lengths, length(args))))
}


##########
##
## Internal function to get variance of binary variable / rule (faster than function sd):
##
var_bin <- function(x) {
  p <- mean(x)
  p*(1L-p)
}


##########
##
## Internal function to check for near equality:
##
is_almost_eq <- function(x, y, tolerance = sqrt(.Machine$double.eps)) {
  stopifnot(is.numeric(x), length(x) == 1L)
  x_abs <- abs(x)
  xy <- if (x_abs > tolerance) {abs(x - y) / x_abs} else {abs(x - y)}
  xy <= tolerance
}


##########
##
## Internal function for transforming tree into a set of rules:
##
# Taken and modified from package partykit, written by Achim Zeileis and 
# Torsten Hothorn
# It has been changed to return all the rules at each node and not just the 
# rules at the terminal nodes. This is done to get the rules as in: 
#   Friedman, J. H., & Popescu, B. E. (2008). Predictive learning via rule 
#   ensembles. The Annals of Applied Statistics, 916-954.


list.rules <- function (x, i = NULL, removecomplements = TRUE, 
                        singleconditions = FALSE, ...) {
  
  ## singleconditions can take values FALSE (rules code only node membership);
  ##   TRUE (rules code both node membership and all individual splits; 
  ##   "only" (rules code only individual splits)
  
  if (is.null(i)) 
    i <- partykit::nodeids(x, terminal = TRUE)
  if (length(i) > 1) {
    # ret <- sapply(i, list.rules, x = x)
    # TODO: Benjamin Christoffersen changed this part. This can be done smarter
    # than finding all and then removing duplicates. I guess the computational
    # cost is low, though
    if (isTRUE(singleconditions)) {
      ret <- lapply(i, list.rules, x = x, simplify = FALSE)
      ret <- c(ret, lapply(i, list.rules, x = x, simplify = FALSE, 
                           singleconditions = TRUE))
    } else if (singleconditions == "only") {
      ret <- lapply(i, list.rules, x = x, simplify = FALSE, 
                           singleconditions = TRUE)
    } else {
      ret <- lapply(i, list.rules, x = x, simplify = FALSE)
    }
    
    ## Find the first rules. We will only keep one of these:
    
    ## TODO: If we apply non-negativity constraints,
    ## the rule that is kept should correlate positively
    ## with the outcome, if we apply negativity constraint,
    ## the rule that is kept should correlate negatively with 
    ## the response.
    ## I.e., if  'lower.limits = 0' or 'upper.limits = 0' was used in calling pre()
    ##
    ## Easier solution may be to just not remove first rule here
    ## E.g., employ rm.firstrule argument (which is true by default)
    if (removecomplements) {
      first_rules <- unique(sapply(ret, "[[", 1))
      first_rule_remove <- first_rules[2]
    }
    
    # Make list of final rules
    ret <- unlist(ret)
    ret <- ret[!duplicated(ret)]
    if (removecomplements) {
      ret <- ret[ret != first_rule_remove]
    }
    ## of rules with only single conditions, retain one of each pair
    if (isTRUE(singleconditions) && removecomplements) {
      mcrs <- ret[grep(" & ", ret)] ## multi-condition rules
      scrs <- ret[-grep(" & ", ret)] ## single-condition rules
      ## keep only odd numbered of those 
      ret <- c(mcrs, scrs[1:length(scrs) %% 2 == 1])
    } else if (singleconditions == "only" && removecomplements) {
      ## keep only odd numbered rules 
      ret <- c(ret[1:length(ret) %% 2 == 1])      
    }
    
    # TODO: this still leaves us with complements for non-terminal rules
    # names(ret) <- if (is.character(i)) 
    #   i else names(x)[i]
    return(ret) # Root node returns here
  }
  
  # Non-root nodes starts here
  if (is.character(i) && !is.null(names(x))) 
    i <- which(names(x) %in% i)
  #stopifnot(length(i) == 1 & is.numeric(i))
  #stopifnot(i <= length(x) & i >= 1)
  i <- as.integer(i)
  # dat <- partykit::data_party(x, i)
  # if (!is.null(x$fitted)) {
  #   findx <- which("(fitted)" == names(dat))[1]
  #   fit <- dat[, findx:ncol(dat), drop = FALSE]
  #   dat <- dat[, -(findx:ncol(dat)), drop = FALSE]
  #   if (ncol(dat) == 0) 
  #     dat <- x$data
  # }
  # else {
  #   fit <- NULL
  #   dat <- x$data
  # }
  dat <- x$data
  rule <- c()
  recFun <- function(node) {
    # if (partykit::id_node(node) == i) {
    #   return(NULL)
    # }
    if (node$id == i) {
      return(NULL)
    }
    # kid <- sapply(partykit::kids_node(node), partykit::id_node)
    kid <- sapply(node$kids, function(x) x$id)
    whichkid <- max(which(kid <= i))
    #split <- partykit::split_node(node)
    split <- node$split
    # ivar <- partykit::varid_split(split)
    ivar <- split$varid
    svar <- names(dat)[ivar]
    # index <- partykit::index_split(split)
    index <- split$index
    if (is.factor(dat[, svar])) {
      # if (is.null(index)) 
      #   index <- ((1:nlevels(dat[, svar])) > partykit::breaks_split(split)) + 1
      if (is.null(index))
        index <- ((1:nlevels(dat[, svar])) > split$breaks) + 1
      slevels <- levels(dat[, svar])[index == whichkid]
      # factor levels not occurring in the node will be coded as NA
      # and should be removed from rule description:
      slevels <- slevels[!is.na(slevels)]
      srule <- paste(svar, " %in% c(\"", paste(slevels, 
                                               collapse = "\", \"", sep = ""), "\")", sep = "")
    } else {
      if (is.null(index)) {
        index <- 1:length(kid)
      }
      # breaks <- cbind(c(-Inf, partykit::breaks_split(split)), c(partykit::breaks_split(split), 
      #                                                           Inf))
      breaks <- cbind(c(-Inf, split$breaks), c(split$breaks, Inf))
      sbreak <- breaks[index == whichkid, ]
      # right <- partykit::right_split(split)
      right <- split$right
      srule <- c()
      if (is.finite(sbreak[1])) {
        srule <- c(srule, paste(svar, ifelse(right, ">", 
                                             ">="), sbreak[1]))
      }
      if (is.finite(sbreak[2])) { 
        srule <- c(srule, paste(svar, ifelse(right, "<=", 
                                             "<"), sbreak[2]))
      }
      srule <- paste(srule, collapse = " & ")
    }
    rule <<- c(rule, srule)
    return(recFun(node[[whichkid]]))
  }
  # node <- recFun(partykit::node_party(x))
  node <- recFun(x$node)
  # paste(rule, collapse = " & ")
  
  if(is.null(rule))
    return(character())
  
  if (isTRUE(singleconditions) || singleconditions == "only") {
    ## keep conditions separate
    rule
  } else {
    ## combine conditions into single rule
    sapply(seq_along(rule), function(r) paste(rule[1:r], collapse = " & "))
  }
  
}





##########
##
## Internal functions for gradient boosting
##
get_intercept_logistic <- function(y, ws = NULL) {
  # # page 484 of:
  # # Bühlmann, Peter, and Torsten Hothorn. "Boosting algorithms: Regularization, 
  # # prediction and model fitting." Statistical Science (2007): 477-505.
  # # or check do the math an figure out that:
  # n <- 1000
  # y <- runif(n) > 1/(1 + exp(-1))
  # w <- runif(n, 0, 2)
  # 
  # glm.fit(
  #   matrix(rep(1, n), ncol = 1), 
  #   y,
  #   family = binomial(), 
  #   weights = w)$coefficients
  # 
  # p <- weighted.mean(y, w)
  # log(p / (1 - p))
  
  p_bar <- if(is.null(ws)) mean(y) else weighted.mean(y, ws)
  log(p_bar / (1 - p_bar))
}


##########
##
## Internal functions for gradient boosting
##
get_y_learn_logistic <- function(eta, y) {
  # See LogitBoost on page 351 of:
  # Friedman, J., Hastie, T., & Tibshirani, R. (2000). Additive logistic 
  # regression: a statistical view of boosting (with discussion and a rejoinder 
  # by the authors). The annals of statistics, 28(2), 337-407.
  
  trunc_fac <- 12
  eta <- pmin(pmax(eta, -trunc_fac), trunc_fac)
  p <- 1 / (1 + exp(-eta))
  (y - p) / sqrt(p * (1 - p))
}


##########
##
## Internal functions for gradient boosting
##
get_intercept_count <- function(y, ws = NULL) {
  lambda_bar <- if(is.null(ws)) mean(y) else weighted.mean(y, ws)
  log(lambda_bar)
}


##########
##
## Internal functions for gradient boosting
##
get_y_learn_count <- function(eta, y) {
  lambda <- exp(eta)
  y - lambda
}


##########
##
## Internal functions for gradient boosting
##
get_intercept_multinomial <- function(y, ws = NULL) {
  p_bar <- if (is.null(ws)) colMeans(y) else apply(y, 2, mean, weights = ws)
  log(p_bar / (1 - p_bar))
}


##########
##
## Internal functions for gradient boosting
##
get_y_learn_multinomial <- function(eta, y) {
  #eta <- cbind(-14:15, 14:-15)
  #y <- cbind(rep(0:1, each = 15), rep(1:0, each = 15))
  eta <- apply(eta, 2, function(eta, trunc_fac = 12) {
    pmin(pmax(eta, -trunc_fac), trunc_fac)
  })
  p <- 1 / (1 + exp(-eta))
  (y - p) / sqrt(p * (1 - p))
}



.get_most_sparse_rule <- function(rules, data){
  #####
  # we could do this faster by evaluating all the rules at once but we do not
  # to reduce the memory usage
  n <- nrow(data)
  sapply(rules, function(r){
    x <- eval(parse(text = r), data)
    if(!is.logical(x))
      stop("non-rulle is passed")
    if(sum(x) > n / 2)
      return(paste0("!(", r, ")"))
    r
  })
}




.get_rules_mat_dense <- function(data, rules){
  if(length(rules) == 0)
    return(NULL)
  
  expr <- parse(text = paste0("cbind(", paste0(rules, collapse = ", "), ")"))
  x <- eval(expr, data)
  colnames(x) <- names(rules)
  x
}

.get_rules_mat_sparse <- function(data, rules){
  if(length(rules) == 0)
    return(NULL)
  
  # See https://stackoverflow.com/a/8844057/5861244. 
  #
  # if all rules where binary then we could use the `lsparseMatrix-classes`. 
  # However, this will not work when we call `glmnet` as it requires a double
  # matrix`
  expr <- paste0("cbind_sparse_vec(", paste0(
    'as(as.numeric(', rules, '), "sparseVector")', collapse = ", "), ")")
  x <- eval(parse(text = expr), data)
  colnames(x) <- names(rules)
  x
}






#' Sampling function generator for specifying varying maximum tree depth 
#' in a prediction rule ensemble (pre)
#' 
#' \code{maxdepth_sampler} generates a random sampling function, governed
#' by a pre-specified average tree depth.
#' 
#' @param av.no.term.nodes integer of length one. Specifies the average 
#' number of terminal nodes in trees used for rule inducation.
#' @param av.tree.depth integer of length one. Specifies the average maximum
#' tree depth in trees used for rule induction.
#' @return Returns a random sampling function with single argument \code{ntrees},
#' which can be supplied to the \code{maxdepth} argument of function 
#' \code{\link{pre}} to specify varying tree depths.
#' @details The original RuleFit implementation varying tree sizes for
#' rule induction. Furthermore, it defined tree size in terms of the number
#' of terminal nodes. In contrast, function \code{\link{pre}} defines the 
#' maximum tree size in terms of a (constant) tree depth. Function 
#' \code{maxdepth_sampler} allows for mimicing the behavior of the
#' orignal RuleFit implementation. In effect, the maximum tree depth is 
#' sampled from an exponential distribution with learning rate 
#' \eqn{1/(\bar{L}-2)}, where \eqn{\bar{L} \ge 2} represents the
#' average number of terminal nodes for trees in the ensemble. See
#' Friedman & Popescu (2008, section 3.3). 
#' @references Friedman, J. H., & Popescu, B. E. (2008). Predictive learning 
#' via rule ensembles. \emph{The Annals of Applied Statistics, 2}(3), 916-954.
#' @seealso \code{\link{pre}}
#' @examples
#' ## RuleFit default is max. 4 terminal nodes, on average:
#' func1 <- maxdepth_sampler()
#' set.seed(42)
#' func1(10)
#' mean(func1(1000))
#' 
#' ## Max. 16 terminal nodes, on average (equals average maxdepth of 4):
#' func2 <- maxdepth_sampler(av.no.term.nodes = 16L)
#' set.seed(42)
#' func2(10)
#' mean(func2(1000))
#' 
#' ## Max. tree depth of 3, on average:
#' func3 <- maxdepth_sampler(av.tree.depth = 3)
#' set.seed(42)
#' func3(10)
#' mean(func3(1000))
#' 
#' ## Max. 2 of terminal nodes, on average (always yields maxdepth of 1):
#' func4 <- maxdepth_sampler(av.no.term.nodes = 2L)
#' set.seed(42)
#' func4(10)
#' mean(func4(1000))
#' 
#' \donttest{## Create rule ensemble with varying maxdepth:
#' set.seed(42)
#' airq.ens <- pre(Ozone ~ ., data = airquality[complete.cases(airquality),],
#'                 maxdepth = func1)
#' airq.ens}
#' @export
maxdepth_sampler <- function(av.no.term.nodes = 4L, av.tree.depth = NULL) {
  function(ntrees, ...) {
    if (!is.null(av.tree.depth)) {
      av.no.term.nodes <- 2^av.tree.depth
    }
    ceiling(log(2 + floor(rexp(ntrees, rate = 1 / (av.no.term.nodes - 2))), base = 2))
  }
}







#' Dealing with rare factor levels in fitting prediction rule ensembles.
#' 
#' Provides a sampling function to be supplied to the \code{sampfrac}
#' argument of function \code{pre}, making sure that each level of specified factor(s)
#' are present in each sample.
#' 
#' @details Categorical predictor variables (factors) with rare levels may be problematic 
#' in boosting algorithms employing sampling (which is employed by default in
#' function \code{pre}).
#' 
#' If a sample in a given boosting iteration does not have any observations with a given
#' (rare) level of a factor, while this level is present in the full training dataset, and 
#' the factor is selected for splitting in the tree, then no prediction for that level of the factor
#' can be generated, resulting in an error. Note that boosting methods other than \code{pre} that also 
#' employ sampling (e.g., \code{gbm} or \code{xgboost}) may not generate an error in such cases, 
#' but also do not document how intermediate predictions are generated in such a case. It is likely that
#' these methods use one-hot-encoding of factors, which from a perspective of model interpretation 
#' introduces new problems, especially when the aim is to obtain a sparse set of rules as in `pre`. 
#'                               
#' With function \code{pre()}, the rare-factor-level issue, if encountered, can be dealt with by the user 
#' in one of the following ways (in random order):
#' 
#' \itemize{
#' \item Use a sampling function that guarantees inclusion of rare factor levels in each sample. E.g., 
#' use \code{rare_level_sampler}, yielding a sampling function which creates training samples 
#' guaranteed to include each level of specified factor(s). Advantage: No loss of information, easy to implement, 
#' guaranteed to solve the issue. Disadvantage: May result in oversampling 
#' of observations with rare factor levels, potentially biasing results. The bias is likely small though, and 
#' will be larger for smaller sample sizes and sampling fractions, and for larger numbers of rare
#' levels. The latter will also increase computational demands. 
#' \item Specify \code{learnrate = 0}. This results in a (su)bagging instead of boosting approach.
#' Advantage: Eliminates the rare-factor-level issue completely, because intermediate predictions
#' need not be computed. Disadvantage: Boosting with low learning rate often improves predictive accuracy.
#' \item Data pre-processing: Before running function \code{pre()}, combine rare factor levels 
#' with other levels of the factors. Advantage: Limited loss of information. Disadvantage: Likely, but 
#' not guaranteed to solve the issue. 
#' \item Data pre-processing: Apply one-hot encoding to the predictor matrix before applying function `pre()`. This can easily be 
#' done through applying function \code{\link[stats]{model.matrix}}. Advantage: Guaranteed to solve the error,
#' easy to implement. Disadvantage: One-hot-encoding increases the number of predictor variables 
#' which may reduce interpretability and, but probably to a lesser extent, accuracy.                                     
#' \item Data pre-processing: Remove observations with rare factor levels from the dataset
#' before running function \code{pre()}. Advantage: Guaranteed to solve the error. Disadvantage: 
#' Removing outliers results in a loss of information, and may bias the results.
#' \item Increase the value of \code{sampfrac} argument of function \code{pre()}. Advantage: Easy to
#' implement. Disadvantage: Larger samples are more likely but not guaranteed to contain all possible 
#' factor levels, thus not guaranteed to solve the issue.
#' }
#' @return A sampling function, which generates sub- or bootstrap samples as usual in function \code{pre}, but 
#' checks if all levels of the specified factor(s) are present and adds observation with those levels if not. 
#' If \code{warning = TRUE}, a warning is issued).
#' 
#' @param factors Character vector with name(s) of factors with rare levels. 
#' @inheritParams pre
#' @param warning logical. Whether a warning should be printed if observations with
#' rare factor levels are added to the training sample of the current iteration.
#' 
#' @examples
#' ## Create dataset with two factors containing rare levels
#' dat <- iris[iris$Species != "versicolor", ]
#' dat <- rbind(dat, iris[iris$Species == "versicolor", ][1:5, ])
#' dat$factor2 <- factor(rep(1:21, times = 5))
#' 
#' ## Set up sampling function
#' samp_func <- rare_level_sampler(c("Species", "factor2"), data = dat, 
#'                                   sampfrac = .51, warning = TRUE)
#' 
#' ## Illustrate what it does                                                                   
#' N <- nrow(dat)
#' wts <- rep(1, times = nrow(dat))
#' set.seed(3)
#' dat[samp_func(n = N, weights = wts), ] # single sample
#' for (i in 1:500) dat[samp_func(n = N, weights = wts), ]
#' warnings() # to illustrate warnings that may occur when fitting a full PRE
#' 
#' ## Illustrate use with function pre:
#' ## (Note: low ntrees value merely to reduce computation time for the example)
#' set.seed(42)
#' # iris.ens <- pre(Petal.Width ~ . , data = dat, ntrees = 20) # would yield error
#' iris.ens <- pre(Petal.Width ~ . , data = dat, ntrees = 20, 
#'   sampfrac = samp_func) # should work
#' @seealso \code{\link{pre}}
#' @export
rare_level_sampler <- function(factors, data, sampfrac = .5, warning = FALSE) {
  if (sampfrac < 1) {
    replace <- FALSE
  }  else if (sampfrac == 1) {
    replace <- TRUE
  }
  ids <- list()
  for (i in factors) {
    ids[[i]] <- lapply(unique(data[ , i]), function(x) which(data[ , i] == x))
    names(ids[[i]]) <- unique(data[ , i])
  }
  ret <- function(n, weights, sampfrac. = sampfrac, replace. = replace, ids. = ids,
                  warning. = warning) {
    ## Generate initial random sample:
    sample_ids <- sample(1:n, size = ceiling(sampfrac.*n), replace = replace., 
                         prob = weights)
    for (i in names(ids.)) {
      ## For every factor, for all levels, check if present in sample
      for (j in names(ids.[[i]])) {
        if (!any(ids[[i]][[j]] %in% sample_ids)) {
          ## If not present, add an observation with that level to the sample
          sample_ids <- c(sample_ids, sample(ids[[i]][[j]], size = 1))
          if (warning.) warning("Added observation with level ", j, " for variable ", i," to current sample.")
        }
      }
    }
    return(sample_ids)
  }
  ret
}
