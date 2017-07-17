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
    
    # Find the first rules. We will only keep one of these
    first_rules <- unique(sapply(ret, "[[", 1))
    first_rule_remove <- first_rules[2]
    
    # Make list of final rules
    ret <- unlist(ret)
    ret <- ret[!duplicated(ret)]
    ret <- ret[ret != first_rule_remove]
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