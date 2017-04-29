# Internal function for transforming tree into a set of rules:
# Taken and modified from package partykit, written by Achim Zeileis and 
# Torsten Hothorn
list.rules <- function (x, i = NULL, ...) 
{
  if (is.null(i)) 
    i <- partykit::nodeids(x, terminal = TRUE)
  if (length(i) > 1) {
    ret <- sapply(i, list.rules, x = x)
    names(ret) <- if (is.character(i)) 
      i
    else names(x)[i]
    return(ret)
  }
  if (is.character(i) && !is.null(names(x))) 
    i <- which(names(x) %in% i)
  #stopifnot(length(i) == 1 & is.numeric(i))
  #stopifnot(i <= length(x) & i >= 1)
  i <- as.integer(i)
  dat <- partykit::data_party(x, i)
  if (!is.null(x$fitted)) {
    findx <- which("(fitted)" == names(dat))[1]
    fit <- dat[, findx:ncol(dat), drop = FALSE]
    dat <- dat[, -(findx:ncol(dat)), drop = FALSE]
    if (ncol(dat) == 0) 
      dat <- x$data
  }
  else {
    fit <- NULL
    dat <- x$data
  }
  rule <- c()
  recFun <- function(node) {
    if (partykit::id_node(node) == i) {
      return(NULL)
    }
    kid <- sapply(partykit::kids_node(node), partykit::id_node)
    whichkid <- max(which(kid <= i))
    split <- partykit::split_node(node)
    ivar <- partykit::varid_split(split)
    svar <- names(dat)[ivar]
    index <- partykit::index_split(split)
    if (is.factor(dat[, svar])) {
      if (is.null(index)) 
        index <- ((1:nlevels(dat[, svar])) > partykit::breaks_split(split)) + 
          1
      slevels <- levels(dat[, svar])[index == whichkid]
      srule <- paste(svar, " %in% c(\"", paste(slevels, 
                                               collapse = "\", \"", sep = ""), "\")", sep = "")
    }
    else {
      if (is.null(index)) {
        index <- 1:length(kid)
      }
      breaks <- cbind(c(-Inf, partykit::breaks_split(split)), c(partykit::breaks_split(split), 
                                                                Inf))
      sbreak <- breaks[index == whichkid, ]
      right <- partykit::right_split(split)
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
  node <- recFun(partykit::node_party(x))
  paste(rule, collapse = " & ")
}

#####
# Wrappers for ctree which returns the needed output to get rules

ctree_setup <- function(
  formula, data, weights, subset, na.action = na.pass, 
  control = ctree_control(...), ytrafo = NULL, scores = NULL, 
  ...
){
  if (missing(data)) 
    data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action"), 
             names(mf), 0)
  mf <- mf[c(1, m)]
  formula <- Formula::Formula(formula)
  mf$formula <- formula
  mf$drop.unused.levels <- FALSE
  mf$na.action <- na.action
  mf[[1]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  response <- names(Formula::model.part(formula, mf, lhs = 1))
  weights <- model.weights(mf)
  dat <- mf[, colnames(mf) != "(weights)"]
  if (!is.null(scores)) {
    for (n in names(scores)) {
      sc <- scores[[n]]
      if (is.ordered(dat[[n]]) && nlevels(dat[[n]]) == 
          length(sc)) {
        attr(dat[[n]], "scores") <- as.numeric(sc)
      }
      else {
        warning("scores for variable ", sQuote(n), " ignored")
      }
    }
  }
  if (is.null(weights)) 
    weights <- rep(1, nrow(mf))
  storage.mode(weights) <- "integer"
  nvar <- sum(!(colnames(dat) %in% response))
  
  control$cfun <- eval(bquote(function(...) {
    if (.(control$teststat == "quad")) 
      p <- .pX2(..., pval = .(control$testtype != "Teststatistic"))
    if (.(control$teststat == "max"))
      p <- .pmaxT(..., pval = .(control$testtype != "Teststatistic"))
    names(p) <- c("statistic", "p.value")
    if (.(control$testtype == "Bonferroni")) 
      p["p.value"] <- p["p.value"] * .(min(nvar, control$mtry))
    crit <- p["statistic"]
    if (.(control$testtype != "Teststatistic")) 
      crit <- p["p.value"]
    c(crit, p)
  }))
  environment(control$cfun) <- environment(ctree)
  
  list(dat = dat, response = response, weights = weights, 
       control = control, ytrafo = ytrafo)
}

ctree_minmal <- function (
  dat, response, weights, control, ytrafo) 
{
  .ctree_fit <- with(environment(ctree), .ctree_fit)
  
  tree <- .ctree_fit(dat, response, weights = weights, ctrl = control, 
                     ytrafo = ytrafo)
  fitted <- data.frame(`(fitted)` = fitted_node(tree, dat), 
                       `(weights)` = weights, check.names = FALSE)
  fitted[[3]] <- dat[, response, drop = length(response) == 1]
  names(fitted)[3] <- "(response)"
  ret <- party(tree, data = dat, fitted = fitted)
               # , info = list(call = match.call(), control = control))
  class(ret) <- c("constparty", class(ret))
  # ret$terms <- terms(mf)
  return(ret)
}