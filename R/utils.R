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
  colnames(rulevars) <- names(rules) <- paste0("rule", 1:length(rules))
  
  ## Remove duplicate rules
  if (removeduplicates) {
    # Remove rules with identical support
    duplicates <- duplicated(rulevars, MARGIN = 2)
    duplicates.removed <- rules[duplicates]
    rulevars <- rulevars[, !duplicates, drop = FALSE]
    rules <- rules[!duplicates]
    
  } else 
    duplicates.removed <- NULL
  
  ## Remove complement rules
  if (removecomplements) {
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
  
  rulevars = if (keep_rulevars && length(rules) > 0)
    rulevars[, names(rules)] else NULL
  
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
duplicated.dgCMatrix <- function (dgCMat, MARGIN) {
  MARGIN <- as.integer(MARGIN)
  n <- nrow(dgCMat)
  p <- ncol(dgCMat)
  J <- rep(1:p, diff(dgCMat@p))
  I <- dgCMat@i + 1
  x <- dgCMat@x
  if (MARGIN == 1L) {
    ## check duplicated rows
    names(x) <- J
    RowLst <- split(x, I)
    is_empty <- setdiff(1:n, I)
    result <- duplicated.default(RowLst)
  } else if (MARGIN == 2L) {
    ## check duplicated columns
    names(x) <- I
    ColLst <- split(x, J)
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


list.rules <- function (x, i = NULL, removecomplements = TRUE, ...) {
  if (is.null(i)) 
    i <- partykit::nodeids(x, terminal = TRUE)
  if (length(i) > 1) {
    # ret <- sapply(i, list.rules, x = x)
    # TODO: Benjamin Christoffersen changed this part. This can be done smarter
    # than finding all and then removing duplicates. I guess the computational
    # cost is low, though
    ret <- lapply(i, list.rules, x = x, simplify = FALSE)
    
    # Find the first rules. We will only keep one of these
    
    ## TODO: If we apply non-negativity constraints,
    ## the rule that is kept should correlate positively
    ## with the outcome, if we apply negativity constraint,
    ## the rule that is kep should correlate negatively with 
    ## the response.
    ## I.e., if  'lower.limits = 0' or 'upper.limits = 0' was used in calling pre()
    ##
    ## Easier solution may be to just not remove first rule here
    ## E..g, employ rm.firstrule argument (which is true by default)
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
    }
    else {
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
  
  sapply(seq_along(rule), function(r) paste(rule[1:r], collapse = " & "))
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