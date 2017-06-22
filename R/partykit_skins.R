# Internal function for transforming tree into a set of rules:
# Taken and modified from package partykit, written by Achim Zeileis and 
# Torsten Hothorn
# It has been changed to return all the rules at each node and not just the 
# rules at the terminal nodes. This is done to get the rules as in: 
#   Friedman, J. H., & Popescu, B. E. (2008). Predictive learning via rule 
#   ensembles. The Annals of Applied Statistics, 916-954.
list.rules <- function (x, i = NULL, ...){
  if (is.null(i)) 
    i <- partykit::nodeids(x, terminal = TRUE)
  if (length(i) > 1) {
    # ret <- sapply(i, list.rules, x = x)
    # TODO: Benjamin Christoffersen changed this part. This can be done smarter
    # then finding all and then removing duplicates. I guess the computational
    # cost is low, though
    ret <- lapply(i, list.rules, x = x, simplify = FALSE)
    ret <- unlist(ret)
    ret <- ret[!duplicated(ret)]
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
  sapply(seq_along(rule), function(r) paste(rule[1:r], collapse = " & "))
}


# This an adaptation of the above that gives rules at each node as in:
#   Friedman, J. H., & Popescu, B. E. (2008). Predictive learning via rule 
#   ensembles. The Annals of Applied Statistics, 916-954
# but without complements
list.all_rules_wo_complements <- function(tree, me, dat, leaves){
  #####
  # Part used for the tree
  if(missing(me)){
    leaves <- partykit::nodeids(tree, terminal = TRUE)
    
    return(list.all_rules_wo_complements(
      tree = tree, me = tree$node, dat = tree$data, leaves = leaves))
  }
  
  #####
  # Part used for nodes
  if (me$id %in% leaves) # we reach leaf
    return(NULL)
  
  kids <- sapply(me$kids, function(x) x$id)
  
  # Get rules from child nodes
  if(length(kids) > 0)
    kid_rules <- lapply(
      partykit::kids_node(me), 
      list.all_rules_wo_complements, 
      tree = tree, dat = dat, leaves = leaves)
  
  # Get rules from this node
  split <- me$split
  ivar <- split$varid
  svar <- names(dat)[ivar]
  index <- split$index
  
  if (is.factor(dat[, svar])) {
    # TODO: test with factor
    if (is.null(index))
      index <- ((1:nlevels(dat[, svar])) > split$breaks) + 1
    slevels <- 
      sapply(sort(unique(index)),function(idx)
        levels(dat[, svar])[index == idx], simplify = FALSE)
    
    my_rules <- sapply(
      slevels, function(lvls) paste(
        svar, " %in% c(\"", paste(
          lvls, collapse = "\", \"", sep = ""), "\")", sep = ""))
    
  } else {
    my_break <- split$breaks
    right <- split$right
    
    my_rules <- vector("character", 2)
    my_rules[1] <- paste(svar, ifelse(right, "<=", "<"), my_break)
    my_rules[2] <- paste(svar, ifelse(right, ">", ">="), my_break)
  }
  
  is_any_kid_a_leaf <- any(sapply(kid_rules, is.null))
  
  my_rules <- sapply(1:2, function(i){
    if(!is.null(kid_rules[[i]])){
      out <-  paste(
        my_rules[i], kid_rules[[i]], sep = " & ")
    } else
      out <- NULL
    
    # We only return the first rule
    if(i == 2)
      return(out)
    
    c(my_rules[i], out)
  }, simplify = FALSE) 
  
  do.call(c, my_rules)
}