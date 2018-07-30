##########
##
## Internal function to remove duplicate and complement rules:
##
delete_duplicates_complements <- function(rules, data, removecomplements = TRUE, 
                                          removeduplicates = TRUE, 
                                          return.dupl.compl = FALSE) {
  
  
  ## Generate rule variables:
  expr <- parse(text = paste0("cbind(", paste0(rules, collapse = ", "), ")"))
  rulevars <- eval(expr, data)
  colnames(rulevars) <- names(rules) <- paste0("rule", 1:length(rules))
  
  ## Remove duplicate rules:
  if (removeduplicates) {
    # Remove rules with identical support:
    duplicates <- duplicated(rulevars, MARGIN = 2)
    duplicates.removed <- rules[duplicates]
    rulevars <- rulevars[, !duplicates, drop = FALSE]
    rules <- rules[!duplicates]
  } else {
    duplicates.removed <- NULL
  }
  
  
  ## Remove complement rules:
  if (removecomplements) {
    # Get variance of each rule:
    vars <- apply(rulevars, 2, var_bin)
    vars_distinct <- sapply(
      unique(vars), function(x) c(x, sum(is_almost_eq(x, vars))))
    complements <- vector(mode = "logical", length(vars))
    for (i in 1:ncol(vars_distinct)) {
      if (vars_distinct[2, i] < 2)
        next
    
      indices <- which(is_almost_eq(vars_distinct[1, i], vars))
    
      for (j in 2:length(indices)) {
        indices_prev <- indices[1:(j - 1)] 
        complements[indices_prev] <- 
          complements[indices_prev] | apply(
            rulevars[, indices_prev, drop = FALSE] != rulevars[, indices[j]], 2, all)
      }
    }  
    complements <- which(complements)
    complements.removed <- rules[complements]
    if (length(complements) > 0) {
      rulevars <- rulevars[, -complements, drop = FALSE]
      rules <- rules[-complements]
    }
  } else {
    complements.removed <- NULL
  }
  
  ## Return results:
  if (return.dupl.compl) {
    return(list(rules = rules, 
                duplicates.removed = duplicates.removed,
                complements.removed = complements.removed))
  } else {
    return(rules)
  }
  
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