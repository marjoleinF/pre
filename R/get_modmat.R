
get_modmat <- function(
  # Pass these if you already have an object
  wins_points = NULL, x_scales = NULL,
  # These should be passed in all calls
  formula, data, rules, type, x_names, winsfrac, normalize, 
  # Response variable names is optional:
  y_names = NULL, sparse = FALSE, rulevars = NULL,
  confirmatory = NULL) {
  
  ## TODO: is next line necessary? Only for converting ordered categorical vars to linear terms?
  ## Is only used twice, in get_rulemat function below.
  ## Perhaps should be conditional on ordinal argument of pre()?
  data_org <- data # needed to evaluate rules later
  
  # convert ordered categorical predictor variables to linear terms:
  data[ , sapply(data, is.ordered)] <- 
    as.numeric(data[ , sapply(data, is.ordered)])
  
  ## Add confirmatory terms:
  if (!is.null(confirmatory)) {
    ## Get confirmatory rules and variables:
    conf_vars <- confirmatory[confirmatory %in% x_names]
    conf_rules <- confirmatory[!(confirmatory %in% conf_vars)]
    if (length(conf_rules) > 0L) {
      ## Evaluate rules (add to model matrix later):
      eval_rules <- function(data, rule_descr) {
        1L * with(data, eval(parse(text = rule_descr)))
      }
      conf_rules <- mapply(FUN = eval_rules, rule_descr = conf_rules, 
                           MoreArgs = list(data = data))
    }
  }
  
  #####
  # Perform winsorizing and normalizing
  ## TODO: Allow for not supplying all variables, but only variables with non-zero importances:
  if (type != "rules" && any(sapply(data[ , x_names], is.numeric))) {
    #####
    # if type is not rules, linear terms should be prepared:
    
    # Winsorize numeric variables (section 5 of F&P(2008)):
    if (winsfrac > 0) {
      
      if (miss_wins_points <- is.null(wins_points)) {
        wins_points <- data.frame(varname = x_names, value = NA, lb = NA, ub = NA, 
                                  stringsAsFactors = FALSE)
      }
      
      j <- 0L
      tol <- sqrt(.Machine$double.eps)
      for (i in x_names) {
        j <- j + 1L
        if (is.numeric(data[[i]])) {
          if (miss_wins_points) {
            lim <- quantile(data[ , i], probs = c(winsfrac, 1 - winsfrac))
            wins_points$value[j] <- paste(lim[1L], "<=", i, "<=", lim[2L])
            lb <- lim[1L]
            ub <- lim[2L]
            if (ub - lb < tol) {
              ## If lower and upper bound are equal, do not winsorize and issue warning:
              warning("Variable ", x_names[j], " will be winsorized employing winsfrac = 0, to prevent reducing the variance of its linear term to 0.", immediate. = TRUE)
              wins_points$lb[j] <- min(data[ , i])
              wins_points$ub[j] <- max(data[ , i])
            } else {
              wins_points$lb[j] <- lb
              wins_points$ub[j] <- ub
              data[ , i][data[ , i] < lb] <- lb
              data[ , i][data[ , i] > ub] <- ub
            }
          } else {
            data[ , i][data[ , i] < wins_points$lb[j]] <- wins_points$lb[j]
            data[ , i][data[ , i] > wins_points$ub[j]] <- wins_points$ub[j]
          }
        }
      }
    }
    
    # normalize numeric variables:
    if (normalize) { 
      # Normalize linear terms (section 5 of F&P08), if there are any:
      needs_scaling <- x_names[sapply(data[ , x_names, drop = FALSE], is.numeric)]
      if (length(needs_scaling) > 0) {
        if (is.null(x_scales)) {
          x_scales <- apply(
            data[ , needs_scaling, drop = FALSE], 2L, sd, na.rm = TRUE) / 0.4
        }
        ## check if variables have zero variance (if so, do not scale):
        tol <- sqrt(.Machine$double.eps)
        almost_zero_var_inds <- which(x_scales < tol)
        if (length(almost_zero_var_inds) > 0) {
          # print warning and set all those x_scales to 1
          warning("A total of ", length(almost_zero_var_inds), " variable(s) (e.g.,", paste0(head(needs_scaling[almost_zero_var_inds]), collapse = ", "), ") have sd < ", tol, " and will not be normalized.")  
          # omit from needs_scaling:
          x_scales[almost_zero_var_inds] <- 1
        }
        data[ , needs_scaling] <- scale(
          data[ , needs_scaling, drop = FALSE], center = FALSE, scale = x_scales)
      }
    }
  }

  ## Combine rules and variables:
  if (length(rules) == 0L) rules <- NULL
  if (type == "linear" || is.null(rules)) {
    if (sparse) {
      mf <- model.frame(Formula(formula), data)
      x <- model.Matrix(terms(mf), mf, sparse = TRUE)
      rm(mf)
      
    } else {
      x <- model.matrix(Formula(formula), data = data)
    }
  } else if (type %in% c("both", "rules") && !is.null(rules)) {
    if (sparse) {
      x <- if (is.null(rulevars)) { 
        .get_rules_mat_sparse(data_org, rules) 
        } else {
          rulevars
        }
      if (type == "both") {
        mf <- model.frame(Formula(formula), data)
        x <- cbind(model.Matrix(terms(mf), mf, sparse = TRUE), x)
        rm(mf)
      }
    } else { 
      x <- if (is.null(rulevars)) {
        .get_rules_mat_dense(data_org, rules)
        } else { 
          rulevars
        }
      if (type == "both") {
        x <- cbind(model.matrix(Formula(formula), data = data), x)
      }
    }
    
  } else { 
    stop("not implemented with type ", sQuote(type), " and is.null(rules) is ", 
         sQuote(is.null(rules)))
  }
  
  #####
  # Remove intercept
  x <- x[, colnames(x) != "(Intercept)", drop = FALSE]

  ## Add confirmatory terms:
  if (!is.null(confirmatory)) {
    ## Get confirmatory rules and variables:
    #conf_vars <- confirmatory[confirmatory %in% x_names]
    #conf_rules <- confirmatory[!(confirmatory %in% conf_vars)]
    if (length(conf_vars) > 0L) {
      ## If type == "rules", add variables to model matrix:
      if (type == "rules") {
        x <- cbind(x, as.matrix(data[conf_vars], rownames.force = FALSE))
        ## TODO: Then check if it should be winsorized or normalized, or issue warning if not
      }
    }
    if (is.matrix(conf_rules)) {
      #if (!sparse) {
        x <- cbind(x, conf_rules)
        ## I believe both sparse and non-sparse matrices can use cbind:
        ## Matrix::cBind returns:
        ## .Defunct(msg = "'cBind' is defunct; 'base::cbind' handles S4 objects since R 3.2.0")
      #} else {
        ## TODO: implement something for sparse matrices.
      #}
    }
  }
  
  if (is.null(y_names)) {
    y <- NULL
  } else {
    y <- data[ , y_names]
    if (is.Surv(y) || length(y_names) > 1L) {
      y <- as.matrix(y)
    }
  }
  
  if (!exists("wins_points", inherits = FALSE)) { wins_points <- NULL }
  if (!exists("x_scales", inherits = FALSE)) { x_scales <- NULL }
  
  list(x = x, y = y, x_scales = x_scales, wins_points = wins_points)
}