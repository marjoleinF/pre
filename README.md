pre: An R package for deriving prediction rule ensembles
========================================================

pre is an R package for deriving prediction rule ensembles for binary and continuous outcome variables. Input variables may be numeric, ordinal and nominal. The package implements the algorithm for deriving prediction rule ensembles as described in Jerome H. Friedman & Popescu (2008), with some improvements and adjustments. The most important improvements and adjustments are:

1.  The package is completely R based, allowing users better accessible results and more control over the parameters used for generating the prediction rule ensemble
2.  The unbiased tree induction algorithm of Hothorn, Hornik, & Zeileis (2006) is used for deriving prediction rules, instead of the classification and regression tree (CART) algorithm, which suffers from biased variable selection.
3.  The package allows for plotting the final rule ensemble as a collection of simple decision trees.
4.  The initial ensemble of prediction rules can be generated as a bagged, boosted and/or random forest ensemble.
5.  Hinge functions of predictor variables may be included as baselearners in the ensemble, as in the multivariate adaptive regression splines technique of Jerome H Friedman (1991).

The pre package is developed to provide useRs a completely R based implementation of the algorithm described by Jerome H. Friedman & Popescu (2008). However, note that pre is under development, and much work still needs to be done. See Fokkema, Smits, Kelderman, & Penninx (2015) for an application of the methods.

Examples
========

To get a first impression of how pre works, we will fit a prediction rule ensemble to predict Ozone levels using the airquality dataset:

``` r
library(pre)
complete <- complete.cases(airquality)
set.seed(42)
airq.ens <- pre(identity(Ozone + Temp) ~ ., data = airquality[complete, ], standardize = TRUE)
```

We can print the resulting ensemble:

``` r
print(airq.ens)
#> 
#> Final ensemble with cv error within 1se of minimum: 
#>   lambda =  6.862571
#>   number of terms = 9
#>   mean cv error (se) = 764.909 (129.2821) 
#> 
#>          rule  coefficient                   description
#>   (Intercept)  141.9909641                          <NA>
#>        rule25  -17.3492833                    Wind > 5.7
#>        rule72   12.2670362   Solar.R > 148 & Wind <= 8.6
#>        rule12  -11.5242457       Wind > 6.3 & Month <= 5
#>        rule62   10.9626070  Solar.R > 149 & Wind <= 10.3
#>       rule104   10.7793652     Wind <= 8 & Solar.R > 201
#>         rule4  -10.6405536                    Wind > 6.3
#>        rule74   -8.5041926    Wind > 6.3 & Solar.R <= 78
#>          Wind   -0.4078209         3.85 <= Wind <= 17.05
#>        rule10   -0.2570948   Wind > 8.6 & Solar.R <= 149
```

We can plot a subsample of the rules (and or/linear terms) in the ensemble:

``` r
plot(airq.ens, max.terms.plot = 9, cex = .75, penalty = "lambda.1se")
```

![](inst/README-figures/README-unnamed-chunk-4-1.png)

We can obtain the estimated coefficients of the ensemble:

``` r
head(coef(airq.ens))
#>           rule coefficient                  description
#> 86 (Intercept)   141.99096                         <NA>
#> 16      rule25   -17.34928                   Wind > 5.7
#> 48      rule72    12.26704  Solar.R > 148 & Wind <= 8.6
#> 9       rule12   -11.52425      Wind > 6.3 & Month <= 5
#> 41      rule62    10.96261 Solar.R > 149 & Wind <= 10.3
#> 75     rule104    10.77937    Wind <= 8 & Solar.R > 201
```

We can assess the importance of the predictor variables, and each of the rules and/or linear terms in the ensemble:

``` r
importance(airq.ens)
```

![](inst/README-figures/README-unnamed-chunk-6-1.png)

    #> $varimps
    #>   varname       imp
    #> 1    Wind 21.142565
    #> 2 Solar.R  9.281859
    #> 3   Month  2.346010
    #> 
    #> $baseimps
    #>      rule                  description       imp coefficient        sd
    #> 1  rule72  Solar.R > 148 & Wind <= 8.6 5.6323040  12.2670362 0.4591414
    #> 2  rule62 Solar.R > 149 & Wind <= 10.3 5.4412030  10.9626070 0.4963421
    #> 3  rule25                   Wind > 5.7 5.4116777 -17.3492833 0.3119251
    #> 4  rule12      Wind > 6.3 & Month <= 5 4.6920208 -11.5242457 0.4071434
    #> 5 rule104    Wind <= 8 & Solar.R > 201 4.1617027  10.7793652 0.3860805
    #> 6   rule4                   Wind > 6.3 3.9398902 -10.6405536 0.3702712
    #> 7  rule74   Wind > 6.3 & Solar.R <= 78 3.2177031  -8.5041926 0.3783667
    #> 8    Wind        3.85 <= Wind <= 17.05 0.1631283  -0.4078209 0.4000000
    #> 9  rule10  Wind > 8.6 & Solar.R <= 149 0.1108044  -0.2570948 0.4309865

We can generate predictions for new observations:

``` r
airq.preds <- predict(airq.ens, newdata = airquality[1:4,])
```

We can assess the effect of a given predictor variable on the outcome through a partial dependence plot:

``` r
singleplot(airq.ens, varname = "Temp")
```

![](inst/README-figures/README-unnamed-chunk-8-1.png)

We can assess the expected prediction error of the ensemble, by default calculated using 10-fold cross validation:

``` r
set.seed(43)
airq.cv <- cvpre(airq.ens)
airq.cv$accuracy
#> $MSE
#> [1] 790.7663
#> 
#> $MAE
#> [1] 22.09351
```

More complex prediction ensembles can be derived with the gpe() function. The abbreviation gpe stands for generalized prediction ensembles, which in addition to rules and linear terms may also include hinge functions of the predictor variables Jerome H Friedman (1991). Addition of hinge functions may improve predictive accuracy (but may also reduce interpretability).

References
==========

Fokkema, M., Smits, N., Kelderman, H., & Penninx, B. W. (2015). Connecting clinical and actuarial prediction with rule-based methods. *Psychological Assessment*, *27*(2), 636.

Friedman, J. H. (1991). Multivariate adaptive regression splines. *The Annals of Statistics*, 1–67.

Friedman, J. H., & Popescu, B. E. (2008). Predictive learning via rule ensembles. *The Annals of Applied Statistics*, *2*(3), 916–954. Retrieved from <http://www.jstor.org/stable/30245114>

Hothorn, T., Hornik, K., & Zeileis, A. (2006). Unbiased recursive partitioning: A conditional inference framework. *Journal of Computational and Graphical Statistics*, *15*(3), 651–674.
