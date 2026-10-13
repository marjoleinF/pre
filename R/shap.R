## TODO: Implement SHAP computation if no rules but only linear terms are selected.

#' Compute SHAP values for a prediction rule ensemble
#'
#' This function computes the marginal (interaction) SHAP values for a prediction rule ensemble as fitted with function
#' \code{pre}.
#'
#' @inheritParams coef.pre
#' @param newdata An optional \code{data.frame} containing observations for which SHAP values should be computed.
#' If \code{NULL}, SHAP values will be computed for the model's training data, obtained from \code{object}.
#' @param reference_data An optional \code{data.frame} containing the reference data used to estimate the expectations 
#' in the SHAP values. If \code{NULL}, the model's training data will be used, obtained from \code{object}.
#' @param interactions A logical value indicating whether marginal interaction SHAP values should be computed,
#' on top of the overall values.
#' @param block_size Computation of SHAP values involves computations on large matrices. For computational reasons, 
#' these are split into small submatrices. This integer denotes the number of rows of the submatrix.
#'
#' @return A list with two objects: 
#' 
#' \code{marginal} An $N$ by $p$ \code{matrix} with SHAP values for each observation (rows) and predictor (columns), 
#' with rownames corresponding to \code{rownames(newdata)}.
#'  
#' \code{interactions} A $p$ by $p$ by $N$ \code{array} with SHAP interaction values. Each slice contains a 
#' $p /times p$ symmetric matrix, of which the diagonal entries present the predictor 
#' variables' main effect SHAP values and the off-diagonal entries represent their interaction SHAP values.
#'   
#' All SHAP values are on the scale of the linear predictor.
#' @export
#' @examples \donttest{## Fit pre, then compute SHAP values for the first 10 training observations 
#' airq <- airquality[complete.cases(airquality), ]
#' set.seed(42) 
#' airq.ens <- pre(Ozone ~ ., data = airq)
#' airq.shap <- shap(airq.ens, newdata = airq[1:10, ])}
#' 
#' @details Adapted code from original written by and used with permission from Giorgio Spadaccini.
#' 
#' @author Giorgio Spadaccini
#' 
#' @seealso \code{\link{print.pre}}, \code{\link{plot.pre}}, 
#' \code{\link{coef.pre}}, \code{\link{importance.pre}}, \code{\link{predict.pre}}, 
#' \code{\link{interact}}, \code{\link{cvpre}} 
#' 
#' 
shap <- function(object, newdata = NULL, reference_data = NULL, 
                 penalty.par.val = "lambda.1se",
                 interactions = FALSE, block_size = 5e3, ...) {
  
  ## Argument checks
  if (!inherits(object, "pre")) stop("Argument object should specify an object of class pre, i.e. an object fitted with function pre.")
  if (object$family %in% c("mgaussian", "multinomial")) stop("SHAP value computation is not yet available for multinomial and multivariate outcomes.")
  if (is.null(newdata)) {
    newdata <- object$data[ , object$x_names]    
  } else {
    if (!all(object$x_names %in% names(newdata))) {
      warning("Argument newdata does not contain all variables used to train the original rule ensemble. An error will likely occur.")
    } else {
      newdata <- newdata[ , object$x_names]
    }
    ## check if factors in data have the same levels as factors in object$data
    if (any(factors <- sapply(object$data[ , object$x_names], is.factor))) {
      for (factor_name in object$x_names[which(factors)]) {
        if (!all(levels(object$data[, factor_name]) == levels(newdata[ , factor_name]))) {
          warning("Feature ", factor_name, " has levels ", paste(levels(newdata[, factor_name]), collapse = ", "), 
                  " in newdata, while ", factor_name, " has levels ", 
                  paste(levels(object$data[, factor_name]), collapse = ", "), 
                  " in the training data used to fit the ensemble. Computed SHAP values may be incorrect, make sure the levels of factors in data match those in the data used to fit the rule ensemble, then rerun.")
        }
      }
    }
  }
  if (is.null(reference_data)) {
    reference_data <- object$data[ , object$x_names]
  } else {
    reference_data <- reference_data[ , object$x_names]    
    if (!all(object$x_names %in% names(reference_data))) {
      warning("Argument reference_data does not contain all variables used to train the original rule ensemble. An error might result.")
    }
    if (any(factors <- sapply(object$data[ , object$x_names], is.factor))) {
      for (factor_name in names(object$data)[which(factors)]) {
        if (!all(levels(object$data[, factor_name]) == levels(reference_data[ , factor_name]))) {
          warning("Feature ", factor_name, " has levels ", paste(levels(reference_data[, factor_name]), collapse = ", "), 
                  "in newdata, while ", factor_name, " has levels ", paste(levels(object$data[, factor_name]), collapse = ", "), 
                  " in the training data used to fit the ensemble. Computed SHAP values may be incorrect, make sure the levels of factors in data match those in the data used to fit the rule ensemble, then rerun.")
        }
      }
    }
  }
  
  ## Compute basic quantities
  p <- length(object$x_names)
  n_test <- nrow(newdata)
  
  ## Get model coefs and select only those rules that are retained in final model
  CoeffVec <- coef(object$glmnet.fit, s = penalty.par.val, ...)
  if (rownames(CoeffVec)[1] != "(Intercept)") stop("Function explain.pre assumes that the glmnet fit contains an intercept but coef(object$glmnet.fit) does not return an intercept. Computations halted.")
  sel_terms <- abs(CoeffVec) > 0
  sel_terms <- c(FALSE, sel_terms[-1L]) ## the fitted object contains an intercept which should not be used  
  sel_rules <- grepv("rule", rownames(object$glmnet.fit$glmnet.fit$beta)[sel_terms[-1L]])
  
  ## If any terms are selected, compute SHAP values
  if (all(!sel_terms)) {
    shapleys_mat <- data.frame(matrix(0L, nrow = n_test, ncol = p))
    rownames(shapleys_mat) <- rownames(newdata)
    names(shapleys_mat) <- object$x_names
    shapleys_inter_df <- NULL
    warning("No rules or linear terms were selected in the object with the penalty.par.val specified. All SHAP values are 0, interaction SHAP values are also all 0 and will not be returned.")
  } else {
    if (length(sel_rules) > 0L) {
      rule_objects <- RuleMats(object$rules[sel_rules, ]$description, newdata)
      if (!is.null(reference_data)) {
        reference_rule_objects <- RuleMats(object$rules[sel_rules, ]$description, reference_data)      
      }
    } else {
      stop("No rules were selected in the final ensemble with specified penalty.parameter.value. SHAP computation not yet implemented. Please contact package author if you need this.")
      rule_objects <- reference_rule_objects <- list(Rs = NULL)
    }
    WeightMatrices <- ShapleyMats(data = reference_data, 
                                  data_test = newdata,
                                  Rs = reference_rule_objects$Rs,
                                  Rs_test = rule_objects$Rs,
                                  id_mat = rule_objects$RulePredMat,
                                  interactions = interactions)
  
    ## WeightMatrices$marginal and WeightMatrices$marginal return matrices that are N*p x 
    ## a column for each selected rule and all linear terms in the model matrix (i.e. all 
    ## linear terms are retained and factors are encoded as nlevels-1L dummy variables) 
    CoeffVec <- CoeffVec[sel_terms, ]
    ## Collect all linear terms and add selected rules
    lin_terms <- colnames(object$modmat)[!grepl("rule", colnames(object$modmat))]
    lin_terms <- setNames(rep(0, times = length(lin_terms)), lin_terms)
    lin_terms[names(lin_terms) %in% names(CoeffVec)] <- CoeffVec[names(CoeffVec) %in% names(lin_terms)]
    rule_terms <- CoeffVec[names(CoeffVec)[grepl("rule", names(CoeffVec))]]
    CoeffVec <- c(lin_terms, rule_terms)
    
    ## List predictor values used for denoting shap interaction values later (should not be winsorized)
    x1 <- c(unlist(lapply(newdata,as.numeric), use.names = FALSE))
    x2 <- c(unlist(lapply(newdata, \(x) rep(as.numeric(x), times = p))))
      
    ## Winsorize data (only linear terms evaluated from here, no rules, so can safely winsorize) 
    if (!all(lin_terms == 0)) {
      for (i in object$wins_points$varname) {
        if (!is.na(object$wins_points[object$wins_points$varname == i, "value"])) {
          newdata[newdata[ , i] < object$wins_points[object$wins_points$varname == i, "lb"], i] <- 
            object$wins_points[object$wins_points$varname == i, "lb"]
          newdata[newdata[ , i] > object$wins_points[object$wins_points$varname == i, "ub"], i] <- 
            object$wins_points[object$wins_points$varname == i, "ub"]
        }
      }
      if (!is.null(reference_data)) {
        for (i in object$wins_points$varname) {
          if (!is.na(object$wins_points[object$wins_points$varname == i, "value"])) {
            reference_data[reference_data[ , i] < object$wins_points[object$wins_points$varname == i, "lb"], i] <- 
              object$wins_points[object$wins_points$varname == i, "lb"]
            reference_data[reference_data[ , i] > object$wins_points[object$wins_points$varname == i, "ub"], i] <- 
              object$wins_points[object$wins_points$varname == i, "ub"]
          }
        }
      }
    }

    ## First do main effects (n x p)
    shapleys_mat <- matrix(WeightMatrices$marginal %*% CoeffVec, 
                           nrow = n_test, ncol = p,
                           dimnames = list(rownames(newdata), object$x_names))
    
    shapleys_main <- matrix( ## n x p    
      WeightMatrices$main %*% CoeffVec,
      nrow = n_test,
      ncol = p,
      dimnames = list(rownames(newdata), object$x_names)
    )
  
    ## Now do interactions, if requested
    if (interactions) {
      vals <- rep(NA, nrow(WeightMatrices$interactions))
      for (i in 1:ceiling(nrow(WeightMatrices$interactions)/block_size)) {
        block_i <- (block_size*(i-1)+1):min(block_size*i,nrow(WeightMatrices$interactions))
        vals[block_i] <- WeightMatrices$interactions[block_i, , drop=FALSE] %*% CoeffVec
      }
      shapleys_inter <- array(0, dim = c(n_test, p, p), ## shapviz wants such an array
                              dimnames = list(rownames(newdata), object$x_names, object$x_names)) 
      
      ## Extract all SHAP values for each row of the original data
      val_ids <- rep_len(1:n_test, length.out  = length(vals))
      for (i in 1:n_test) {
        shapleys_inter[i , , ] <- vals[val_ids == i]
        diag(shapleys_inter[i , , ]) <- shapleys_main[i, ]
      }
      
      #shapleys_inter_df <- data.frame(obs_id = rownames(newdata),
      #                                predictor1 = rep(object$x_names, each=n_test), x1 = x1, 
      #                                predictor2 = rep(object$x_names, each=n_test*p), x2 = x2,
      #                                value = c(vals))
      #rownames(shapleys_inter_df) <- NULL
    } else {
      #shapleys_inter_df <- NULL
      shapleys_inter <- NULL
    }
  }
  
  return(list(marginal=shapleys_mat, interactions=shapleys_inter))
}




#' Compute a matrix with SHAP values of each rule
#'
#' This function computes the SHAP values of each rule and
#' arranges them in a matrix that can be used to compute SHAP values of the
#' model as a whole.
#'
#' @param data a dataframe containing the data used to estimate the expectations in SHAP values
#' @param data_test a dataframe containing the points to compute the SHAP values of. By default,
#' this coincides with the datapoints used to estimate SHAP values.
#' @param Rs a list of as many matrices as there are rules to compute the SHAP values for.
#' The \eqn{j}-th element is a matrix corresponding to the \eqn{j}-th rule. Each of its columns
#' corresponds to the 0-1 encoding of a subrule of the \eqn{j}-th rule, as observed in the data
#' provided with input parameter \code{data}. Can be computed with the \link{RuleMats} function.
#' @param Rs_test same as Rs, but computed for the (possibly different) observations provided
#' from the \code{data_test} parameter.
#' @param id_mat A matrix with as many rows as there are rules and as many columns as there
#' are predictors. The \eqn{j}-th column has ones on the entries corresponding to rules that
#' involve the \eqn{j}-th predictor, while all remaining entries are zeroes. Can be computed
#' with the \link{RuleMats} function.
#' @param interactions A logical parameter determining whether interaction SHAP values should
#' also be computed.
#' @return marginal A matrix with \eqn{n\cdot p} rows and as many columns as
#' there are terms (both linear and rules). It is obtained by vertically
#' stacking matrices of \eqn{n} rows. Each submatrix focuses on the SHAP values
#' of a different predictor: the \eqn{(i,k)}-th entry of the \eqn{j}-th of
#' such submatrices represents the contribution of the \eqn{k}-th term to
#' the SHAP value of the \eqn{i}-th datapoint for the \eqn{j}-th predictor.
#' The first \eqn{p} terms are the linear terms, and the remaining columns
#' refer to the rules.
#' @return interaction A matrix with \eqn{n\cdot p^2} rows and as many columns as
#' there are terms (both linear and rules). It is obtained by vertically
#' stacking matrices of \eqn{n \cdot p} rows. Each submatrix is in turn split
#' into \eqn{p} subsubmatrices which focuses on the interaction SHAP values
#' of a different pair of predictors: the \eqn{(i,k)}-th entry of the \eqn{j}-th
#' subsubmatrix of the \eqn{j'}-th submatrix represents the contribution of the
#' \eqn{k}-th term to the SHAP value of the \eqn{i}-th datapoint for the
#' interaction between the \eqn{j}-th and the \eqn{j'}-th predictor.
#' The first \eqn{p} terms are the linear terms, and the remaining columns
#' refer to the rules.
#' 
#' @author Giorgio Spadaccini
#' 
#' @details Code written by and used with permission from Giorgio Spadaccini.
ShapleyMats <- function(data, data_test=data, Rs, Rs_test=Rs, id_mat, interactions=FALSE) {
  ## Obtain basic quantities
  n=nrow(data)
  n_test=nrow(data_test)
  X=model.matrix(~ ., data)[,-1] #remove intercept, model.matrix adds it and doing ~ .-1 does not fix it
  X_test=model.matrix(~ ., data_test)[,-1]
  p=ncol(data)
  pX=ncol(X)
  q=nrow(id_mat)
  P=pX+q
  
  #compute how many columns each predictor takes up
  ndummies=unlist(lapply(data,nlevels))-1
  ndummies[ndummies < 0]=1
  
  #To avoid the computational cost of repeatedly updating sparse Matrices, i.e.
  #their dataframes, store everything as preallocated vectors i,j,x (SparseDataFrame).
  #in marginal/main effects, each rule produces n_test*d non-zero coefficients
  #in interactions case, each rule produces n_test*d*(d-1) non-zero coefficients
  #for each of these, use a pointer to keep track of the first entry available
  #for writing (all i,j,x of the same type share same pointer)
  if (is.null(Rs)) {
    d_vec <- pX
  } else {
    d_vec=sapply(Rs,ncol)
    i_marg=j_marg=x_marg=integer(sum(d_vec)*n_test)
    i_inter=j_inter=x_inter=integer(sum(d_vec*(d_vec-1))*n_test)
    first_wrt_marg=first_wrt_inter=1
  }
  
  
  #The first p columns are for the linear terms:
  i_marg[1:(n_test*pX)]=(rep(1:p,times=ndummies*n_test)-1)*n_test+1:n_test
  j_marg[1:(n_test*pX)]=rep(1:pX,each=n_test)
  x_marg[1:(n_test*pX)]=c(X_test)-rep(colMeans(X),each=n_test)
  first_wrt_marg=n_test*pX+1
  
  
  #The remaining ones are for the rules:
  #create a progress bar
  #pb = utils::txtProgressBar(min = 0, max = length(Rs), initial = 0,style=3)
  
  #Go through all the rule contributions one by one
  for(i in 1:length(Rs)) {
    d=d_vec[i]
    #If d=1, then the contribution is like for linear terms (and interactions=0)
    if(d==1){
      #Involved predictor
      inv_pred=which(id_mat[i,]==1)
      
      #Update matrix
      i_marg[seq(first_wrt_marg,length.out=n_test)]=
        (inv_pred-1)*n_test+1:n_test
      j_marg[seq(first_wrt_marg,length.out=n_test)]=pX+i
      x_marg[seq(first_wrt_marg,length.out=n_test)]=
        c(Rs_test[[i]])-mean(Rs[[i]])
      
      #Update pointer
      first_wrt_marg=first_wrt_marg+n_test
      
      #Update progress bar
      #utils::setTxtProgressBar(pb,i)
      next
    }
    
    #Otherwise, if d>1, we compute the Shapleys properly
    SSt=tcrossprod(!Rs_test[[i]],!Rs[[i]])
    SSt_checks=which(SSt==0,arr.ind = TRUE)
    SSt_points=data.frame(x=SSt_checks[,1],y=SSt_checks[,2])
    
    #Now go through all involved predictors and update the contributions
    inv_pred=which(id_mat[i,]==1)
    for(j in 1:d){
      #Start with updating the marginal shapleys
      
      #check the extra condition R_i(x_i)R_i(y_i)=0. These points are contributing
      contributions=SSt_points[(Rs_test[[i]][SSt_points$x,j]*Rs[[i]][SSt_points$y,j])==0,]
      
      #For these points, compute the weights. Arrange them in a matrix
      Weights=matrix(0,nrow=n_test,ncol=n)
      qx=rowSums(Rs_test[[i]][,-j,drop=F])
      qy=rowSums(Rs[[i]][,-j,drop=F])
      Weights[as.matrix(contributions)]=
        (Rs_test[[i]][contributions$x,j]-Rs[[i]][contributions$y,j])/
        (n*(d-qx[contributions$x])*
           choose(2*d-qx[contributions$x]-qy[contributions$y]-1,d-qx[contributions$x]))
      
      sum_contributions=rowSums(Weights)
      
      #Add the weights to the matrix
      #Update matrix
      i_marg[seq(first_wrt_marg,length.out=n_test)]=
        (inv_pred[j]-1)*n_test+1:n_test
      j_marg[seq(first_wrt_marg,length.out=n_test)]=pX+i
      x_marg[seq(first_wrt_marg,length.out=n_test)]=
        sum_contributions
      #Update pointer
      first_wrt_marg=first_wrt_marg+n_test
      
      
      #Now update the interaction shapleys, if so requested
      if(interactions){
        #(to avoid double counting, only do predictors before j)
        for(k in seq(1,length.out=j-1)){
          #Compute interaction-specific q(x)
          qx_int=qx-Rs_test[[i]][,k]
          qy_int=qy-Rs[[i]][,k]
          
          #Compute interaction specific weights (you can use contributions df, since
          #R_j(x_j)=R_j(y_j) implies no contribution to interactions)
          Weights=matrix(0,nrow=n_test,ncol=n)
          Weights[as.matrix(contributions)]=
            (Rs_test[[i]][contributions$x,j]*Rs_test[[i]][contributions$x,k]
             +Rs[[i]][contributions$y,j]*Rs[[i]][contributions$y,k]
             -Rs[[i]][contributions$y,j]*Rs_test[[i]][contributions$x,k]
             -Rs_test[[i]][contributions$x,j]*Rs[[i]][contributions$y,k])/
            (2*n*(d-1-qx_int[contributions$x])*
               choose(2*d-qx_int[contributions$x]-qy_int[contributions$y]-3,d-1-qx_int[contributions$x]))
          
          sum_contributions=rowSums(Weights)
          
          #Add the weights to the interaction matrix (both as j,k and k,j
          #since we wrote the combination only as k<j)
          #Update matrix
          i_inter[seq(first_wrt_inter,length.out=2*n_test)]=
            c((inv_pred[j]-1)*n_test*p+(inv_pred[k]-1)*n_test+1:n_test,
              (inv_pred[k]-1)*n_test*p+(inv_pred[j]-1)*n_test+1:n_test)
          j_inter[seq(first_wrt_inter,length.out=2*n_test)]=pX+i
          x_inter[seq(first_wrt_inter,length.out=2*n_test)]=
            sum_contributions #this whole vector will be repeated twice
          #Update pointer
          first_wrt_inter=first_wrt_inter+2*n_test
        }
      }
      
    }
    
    #Update progress bar
    #utils::setTxtProgressBar(pb,i)
  }
  
  #close progress bar
  #close(pb)
  
  #Using i,j,x computed above, build the sparse matrices
  #(if contributions$x does not contain all points, then i,j,x are shorter
  #on the off chance that this happens, shorten the vectors to prevent error)
  
  Shapleymat=Matrix::sparseMatrix(i=i_marg[seq(1,length.out=first_wrt_marg-1)],
                                  j=j_marg[seq(1,length.out=first_wrt_marg-1)],
                                  x=x_marg[seq(1,length.out=first_wrt_marg-1)],
                                  dims=c(n_test*p,P))
  
  InterShapleymat=Matrix::sparseMatrix(i=i_inter[seq(1,length.out=first_wrt_inter-1)],
                                       j=j_inter[seq(1,length.out=first_wrt_inter-1)],
                                       x=x_inter[seq(1,length.out=first_wrt_inter-1)],
                                       dims=c(n_test*p^2,P))
  
  #The interaction Shapley matrix is missing the main effects. Compute them by
  #deducting the interactions from the marginal shapleys
  MainShapleymat=InterShapleymat
  dim(MainShapleymat)=c(n_test,p^2*P)
  #Prepare matrix to perform block sum
  blocksum=Matrix::sparseMatrix(i=1:(p*P),j=1:(p*P),x=1) %x% matrix(1,nrow=p,ncol=1)
  #Sum up all interactions per predictor, per point
  MainShapleymat=MainShapleymat %*%blocksum
  #Put them back in the same format as marginal shapley values
  dim(MainShapleymat)=c(n_test*p,P)
  #Subtract the sums of interactions from the marginal shapleys. That's the main effect
  MainShapleymat = Shapleymat-MainShapleymat
  
  #Return the matrix
  return(list(marginal=Shapleymat,main=MainShapleymat,
              interactions=InterShapleymat))
}




#' Compute matrices that identify the rules
#'
#' This function computes the matrices that are needed to compute marginal
#' Shapley values of a RuleSHAP model. These matrices check which points satisfy
#' each of the subrules of each rule. They also check which predictors are involved
#' in the definition of each rule.
#'
#' @param rules a character vector containing the rules to compute the matrices for
#' @param x_df a dataframe containing the observations to compute the matrices for
#'
#' @return RulePredMat A matrix with as many rows as there are rules and as many
#' columns as there are predictors. The \eqn{j}-th column has ones on the entries
#' corresponding to rules that involve the \eqn{j}-th predictor, while all
#' remaining entries are zeroes. Can be computed with the \link{RuleMats} function.
#' @return Rs a list of as many matrices as there are rules. The \eqn{j}-th element
#' is a matrix corresponding to the \eqn{j}-th rule. Each of its columns corresponds
#' to the 0-1 encoding of a subrule of the \eqn{j}-th rule, as observed in the data
#' provided with input parameter \code{x}.
#' 
#' @author Giorgio Spadaccini
RuleMats <- function(rules, x_df) {
  rules_separated=stringr::str_split_fixed(rules," & ",n=1+max(stringr::str_count(rules,'&')))
  vars_in_term <- gsub( " .*$", "",rules_separated)
  
  #Codify rules_separated into coding language
  Rulesmat=matrix(paste0("x_df$",rules_separated),ncol=ncol(rules_separated))
  Rulesmat[Rulesmat=="x_df$"]=""
  
  #Some subrules might involve the same predictor. Those, we'd still like together
  for(i in 1:nrow(vars_in_term)){
    for(x_name in unique(vars_in_term[i,])){
      #unique(vars_in_term[i,]) also includes the name "". Skip that case
      if(x_name==''){
        next
      }
      
      #Find all spots sharing the same predictor x_name
      indices=which(x_name==vars_in_term[i,])
      
      #Join back these subrules in rules_separated, place them in first occurrance
      rules_separated[i,indices[1]]=paste(rules_separated[i,indices],collapse=' & ')
      #In the remaining spots, we need to delete the subrule, it was already merged
      rules_separated[i,indices[-1]]=''
      
      #Do the same with Rulesmat
      Rulesmat[i,indices[1]]=paste(Rulesmat[i,indices],collapse=' & ')
      Rulesmat[i,indices[-1]]=''
      
      #Do the same with vars_in_term (nothing to collapse, only delete copies)
      vars_in_term[i,indices[-1]]=''
    }
  }
  #Define a matrix where M_{i,j} tells if x_j is involved in the i-th rule
  x_names <- names(x_df)
  id_mat <- t(apply(vars_in_term,MARGIN=1,FUN=function(x){x_names %in% x}))
  colnames(id_mat)=x_names
  
  #Rulesmat defined the rules in coding. Create a function that runs it
  parseval=function(text){return(eval(parse(text=text)))}
  
  #Use this function to create a list of matrices. Each matrix is \{R_i(x^{(j)}_i)\}_{i,j}
  Rs=list()
  for(i in 1:nrow(Rulesmat)){
    #Build matrix R for the i-th rule
    R=matrix(unlist(apply(Rulesmat[i,,drop=FALSE],MARGIN=2,FUN = parseval)),
             nrow=nrow(x_df))
    colnames(R)=vars_in_term[i,vars_in_term[i,]!='']
    
    #Using vars_in_mat, reorder the matrices by predictor
    Rs[[i]]=R[,order(match(colnames(R),x_names)),drop=FALSE]
  }
  
  return(list(RulePredMat=id_mat,Rs=Rs))
}
