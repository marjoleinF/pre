# pre version 1.0.8 (2025-..-..)

* Function pre now allows use of function randomForest (from package of 
  same name) for rule induction though argument randomForest.

* Functions singleplot and pairplot now plot observed values for 
  predictor(s) of interest as a rug on the x-axis, and invisibly return 
  the evaluated partial dependence function(s) for further use.




# pre version 1.0.7 (2024-11-01)

* Bug fixed in computation of importances: factors included as linear
  terms did not contribute to variable importances. Now they do.
  
* Function mi_pre added, to allow for fitting a single prediction rule
  ensemble on multiply imputed datasets.

* Function prune_pre added, to extract optimal penalty values for an
  ensemble of given size.
  
* Adaptive lasso has been implemented in function pre. Arguments ad_alpha 
  and ad_penalty allow for fitting the final rule ensemble using adaptive 
  lasso. See documentation of function pre and vignette("relaxed", "pre").

* Partial dependence plotting functions singleplot and pairplot now
  also support multinomial and multivariate outcomes.

* Vignette added detailing how function pre's computatation time can be 
  reduced.
 
 
  

# pre version 1.0.6 (2023-02-13)

* No major changes.




# pre version 1.0.5 (2022-06-10)

* No major changes, several bug fixes, dependence on package akima 
  eliminated.




# pre version 1.0.4 (2022-03-30)

* Implemented relaxed lasso fits and added vignette on relaxed fits, 
  see vignette("relaxed", "pre").

* More informative warning message if no rules could be derived.




# pre version 1.0.3 (2022-02-24)

* No news, update to fix failing test.




# pre version 1.0.2 (2022-02-24)

* Added vignette explaining how to tune parameters of function pre() using
  package caret.

* Added vignette on how to deal with missing values.

* Added support for supplying a tibble instead of a data.frame to function pre().

* Function importance is now an S3method.

* Argument cex.axis of function (method) importance is now passed correctly.

* Argument weights of functions pre() and gpe() now correctly passed to rule 
  induction functions.

* Added progress bar indicating tree fitting progress.

* Argument singleconditions implemented in function pre() (experimental).




# pre version 1.0.1 (2021-07-04)

* Added function rare_level_sampler to deal with errors related to
  new factor levels.




# pre version 1.0.0 (2020-05-03)

* Minor changes for compatibility with changed data.frame function 
  in R version 4.0.0.

* References to paper in Journal of Statistical Software included.




# pre version 0.7.2 (2019-11-11)

* Only minor internal changes, required for compatibility with new
  version of glmnet package.




# pre version 0.7.1 (2019-04-24)

* Bugs fixed in caret_pre_model: tuning of penalty.par.val argument
  always yielded results for "lambda.1se" only. Results are now
  correctly returned for "lambda.min" and "lambda.1se". caret's 
  varImp() and predictors() not supported (perhaps temporarily), 
  as these would always employ default penalty.par.val of
  "lambda.1se". 

* Bugs fixed in explain(). 




# pre version 0.7 (2019-03-30)

* Added support for sparse rule matrix, which can be invoked through sparse 
  argument in pre(). If sparse = TRUE, memory usage will be reduced and 
  computation speed may be improved for large datasets.
 
* Added function explain(), which provides (graphical) explanations of the 
  ensemble's predictions at the individual observation level.




# pre version 0.6 (2018-08-03)

* Added support for survival responses (i.e., family = "cox") in pre() 

* Added summary methods for pre and gpe. 

* Extended support to all response variable types available in pre() for 
  functions plot(), importance() and cvpre().
 
* plot.pre now allows for specifying separate plotting colors for rules 
  with positive and negative coefficients.

* coef and print methods for pre now return descriptions for the intercept 
  (and factor variables), thanks to suggestion by Stephen Milborrow.

* Bug fix in pre(): ordered factors no longer yield error. Implemented 
  new argument 'ordinal' in pre(), which specifies how ordered factors 
  should be processed.

* Bug fix in cvpre(): pclass argument now processed correctly.

* Bug fix in cvpre(): previously, SDs insteas of SEs were returned for 
  binary classification. Accurate standard errors are returned now.

* Bugs fixed in coef.pre(), print.pre(), plot.pre() and importance() 
  when tree.unbiased = FALSE, thanks to a bug report by Stephen 
  Milborrow. 
 



# pre version 0.5 (2018-05-07)

* Function pre() now also supports multinomial and multivariate 
  gaussian response variables.
  
* Function pre() now has argument 'tree.unbiased'; if set to FALSE, 
  the CART algorithm (as implemented in package 'rpart') is employed 
  for rule induction.
  
* Argument 'maxdepth' of function pre() allows for specifying 
  varying maximum depth across trees, through specifying a 
  vector of length ntrees, or a random number generating function.
  See ?maxdepth.sampler for details.




# pre version 0.4 (2017-08-31)

* Added dataset 'carrillo'
  
* By default, a gradient boosting approach is now taken for all response 
  types. That is, partykit::ctree() and a learning rate of .01 is 
  employed by default. Alternatively, glmtree() can be employed for tree 
  induction by sprecifying use.grad = FALSE.
  
* The 'family' argument in pre() now takes character strings as well as glm 
  family objects. 
  
* Functions pairplot() and interact() now use HCL instead of highly saturated
  HSV colors as default plotting colors.
  
* Bug fixed in plot.pre: Node directions are now in accordance with 
  rule definition.

* Bug fixed in predict.pre: No error printed when response variable is not 
  supplied.




# pre version 0.3 (2017-08-03):

* Function gpe() added, which fits general prediction ensembles. By default,
  it fits an ensemble of rules, linear and hinge functions. Function gpe()
  allows for specifying custom baselearner generating functions and a custom
  fitting function for the final model.
  
* Numerous bugs fixed, yielding faster computation times and clearer plots 
  with more customization options. 
  
* Added support for count responses. Function pre() now has a 'family' 
  argument, which should be set to 'poisson' for count outcomes (the
  'family' argument is set automatically to 'gaussian' for numeric response
  variables and to 'binomial' for binary response variables (factors)).

* A gradient boosting approach for binary outcomes is applied, by default,
  substantially reducing computation times. This can be turned off through 
  the 'use.grad' argument in function pre().

* The default of the 'learnrate' argument of function pre() has been changed 
  to .01, by default. Before, it was .01 for continuous outcomes, but 0 for
  binary outcomes, to reduce computation time. With gradient boosting
  implemented, computation time is much reduced.
  
* Argument 'tree.control' in function pre() allows for passing arguments to
  partykit tree fitting functions. 
    
* Arguments for the cv.glmnet() function are directly passed through better 
  use of ellipsis (...). Most importantly, this means that argument 'mod.sel.crit' 
  cannot be used anymore and should be referred to as 'type.measure' 
  (which will be directly passed to cv.glmnet). Similarly, 'thres' and 
  'standardize' are not explicit arguments of function pre() anymore and 
  can now be directly passed to cv.glmnet() using ellipsis (...).
      
* Better use of sample weights: weights specified with the 'weights' argument 
  in pre() are now used as weights in the subsampling procedure, instead of
  as observation weights in the tree-fitting procedure.
      
* Added corplot() function, which shows the correlation between the 
  baselearners in the ensemble.

* Function pairplot() returns a heatmap by default, a 3D or contour plot can
  also be requested.
  
* Appearance of plot resulting from interaction() improved.
  



# pre version 0.2 (2017-04-25):

* Added print() and plot() method for objects of class pre.
  
* Added support for using functions like factor() and log() in formula 
  statement of function pre(). (thanks to Bill Venables for suggesting this)
  
* Added support for parallel computating in functions pre(), cvpre(), 
  bsnullinteract() and interact(). 
  
* Winsorizing points used for the linear terms are reported in the description 
  of the base learners returned by coef() and importance(). (Thanks to
  Rishi Sadhir for suggesting this)
  
* Added README file.
  
* Legend included in plot for interaction test statistics.
  
* Fixed importance() function to allow for selecting final ensemble with 
  different value than 'lambda.1se'.
  
* Cleaned up all occurrences of set.seed()
  
* Fixed cvpre() function: penalty.par.val argument now included
  
* Many minor bug fixes.




# pre version 0.1 (2016-12-23):

* First CRAN release.
