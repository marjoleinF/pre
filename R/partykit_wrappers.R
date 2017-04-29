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
  dat, response, control, ytrafo, ...) 
{
  .ctree_fit <- with(environment(ctree), .ctree_fit)
  
  weights <- rep(1, nrow(dat))
  tree <- .ctree_fit(dat, response, weights = weights, ctrl = control, 
                     ytrafo = ytrafo)
  
  # fitted <- data.frame(`(fitted)` = fitted_node(tree, dat),
  #                      `(weights)` = weights, check.names = FALSE)
  # fitted[[3]] <- dat[, response, drop = length(response) == 1]
  # names(fitted)[3] <- "(response)"
  
  # We compute the outcome in each node to start with to reduce the computation
  # time in predict. I guess this is not done because we lose some information 
  # here. E.g. we cannot get predicted probabilities, densties etc. with 
  # predict
  
  end_nodes <- fitted_node(tree, dat)
  resps <- dat[[response]]
  
  # We assume that we either have factors or numeric is used
  FUN <- if(is.numeric(resps))
    with(environment(ctree), .pred_numeric_response) else
      with(environment(ctree), .pred_factor_response)
  
  # We assume a weight of one
  FUN_wrap <- function(y) FUN(y, rep(1, length(y)))
  fits <- tapply(resps, end_nodes, FUN_wrap)
  
  fitted <- data.frame(
    as.integer(names(fits)),
    rep(1, length(fits)),
    fits)
  names(fitted) <- c("(fitted)", "(weights)", "(response)")
  
  ret <- party_minimal(tree, data = dat, fitted)
  # ret <- party(tree, data = dat, fitted = fitted
               # , info = list(call = match.call(), control = control))
  class(ret) <- c("constparty", class(ret))
  # ret$terms <- terms(mf)
  return(ret)
}

party_minimal <- function (
  node, data, fitted = NULL, ...) {
  # stopifnot(inherits(node, "partynode"))
  # stopifnot(inherits(data, "data.frame"))
  # ids <- nodeids(node)[!nodeids(node) %in% nodeids(node, terminal = TRUE)]
  # varids <- unique(unlist(nodeapply(node, ids = ids, FUN = function(x) varid_split(split_node(x)))))
  # stopifnot(varids %in% 1:ncol(data))
  # if (!is.null(fitted)) {
  #   stopifnot(inherits(fitted, "data.frame"))
  #   stopifnot(nrow(data) == 0L | nrow(data) == nrow(fitted))
  #   if (nrow(data) > 0L) {
  #     if (!("(fitted)" %in% names(fitted))) 
  #       fitted[["(fitted)"]] <- fitted_node(node, data = data)
  #   }
  #   else {
  #     stopifnot("(fitted)" == names(fitted)[1L])
  #   }
  #   nt <- nodeids(node, terminal = TRUE)
  #   stopifnot(all(fitted[["(fitted)"]] %in% nt))
  #   node <- as.partynode(node, from = 1L)
  #   nt2 <- nodeids(node, terminal = TRUE)
  #   fitted[["(fitted)"]] <- nt2[match(fitted[["(fitted)"]], 
  #                                     nt)]
  # }
  # else {
  #   node <- as.partynode(node, from = 1L)
  #   if (nrow(data) > 0L & missing(fitted)) 
  #     fitted <- data.frame(`(fitted)` = fitted_node(node, 
  #                                                   data = data), check.names = FALSE)
  # }
  node <- as.partynode(node, from = 1L)
  party <- list(node = node, data = data, fitted = fitted, 
                terms = NULL, names = NULL, info = NULL)
  class(party) <- "party"
  # if (!is.null(terms)) {
  #   stopifnot(inherits(terms, "terms"))
  #   party$terms <- terms
  # }
  # if (!is.null(names)) {
  #   n <- length(nodeids(party, terminal = FALSE))
  #   if (length(names) != n) 
  #     stop("invalid", " ", sQuote("names"), " ", "argument")
  #   party$names <- names
  # }
  party
}


predict_party_minimal <- function (object, newdata = NULL, perm = NULL, ...) 
{
  fitted <- if (is.null(newdata)) {
    object$fitted[["(fitted)"]]
  }
  else {
    terminal <- nodeids(object, terminal = TRUE)
    if (max(terminal) == 1L) {
      rep.int(1L, NROW(newdata))
    }
    else {
      inner <- 1L:max(terminal)
      inner <- inner[-terminal]
      primary_vars <- nodeapply(object, ids = inner, by_node = TRUE, 
                                FUN = function(node) {
                                  varid_split(split_node(node))
                                })
      surrogate_vars <- nodeapply(object, ids = inner, 
                                  by_node = TRUE, FUN = function(node) {
                                    surr <- surrogates_node(node)
                                    if (is.null(surr)) 
                                      return(NULL)
                                    else return(sapply(surr, varid_split))
                                  })
      vnames <- names(object$data)
      if (!is.null(perm)) {
        stopifnot(all(perm %in% vnames))
        perm <- match(perm, vnames)
      }
      unames <- vnames[unique(unlist(c(primary_vars, surrogate_vars)))]
      vclass <- structure(lapply(object$data, class), 
                          .Names = vnames)
      ndnames <- names(newdata)
      ndclass <- structure(lapply(newdata, class), .Names = ndnames)
      checkclass <- all(sapply(unames, function(x) isTRUE(all.equal(vclass[[x]], 
                                                                    ndclass[[x]]))))
      factors <- sapply(unames, function(x) inherits(object$data[[x]], 
                                                     "factor"))
      checkfactors <- all(sapply(unames[factors], function(x) isTRUE(all.equal(levels(object$data[[x]]), 
                                                                               levels(newdata[[x]])))))
      if (all(unames %in% ndnames) && checkclass && checkfactors) {
        vmatch <- match(vnames, ndnames)
        fitted_node(node_party(object), data = newdata, 
                    vmatch = vmatch, perm = perm)
      }
      else {
        if (!is.null(object$terms)) {
          mf <- model.frame(delete.response(object$terms), 
                            newdata)
          fitted_node(node_party(object), data = mf, 
                      vmatch = match(vnames, names(mf)), perm = perm)
        }
        else stop("")
      }
    }
  }
  # predict_party
  
  # Assume that the fitted element has exactly one match with the right value
  # in reponse for the end_node given in each element of the vector  fitted
  object$fitted[["(response)"]][match(fitted, object$fitted[["(fitted)"]])]
}