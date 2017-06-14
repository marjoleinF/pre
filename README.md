pre: An R package for deriving prediction rule ensembles
========================================================

pre is an R package for deriving prediction rule ensembles for binary and continuous outcome variables. Input variables may be numeric, ordinal and nominal. The package implements the algorithm for deriving prediction rule ensembles as described in Friedman & Popescu (2008), with some improvements and adjustments. The most important improvements and adjustments are:

1.  The pre package is completely R based, allowing users better access to the results and more control over the parameters used for generating the prediction rule ensemble
2.  An unbiased tree induction algorithm is used for deriving prediction rules. Friedman & Popescu (2008) employed the classification and regression tree (CART) algorithm, but this suffers from biased variable selection.
3.  The package allows for plotting the final rule ensemble as a collection of simple decision trees.
4.  The initial ensemble of prediction rules can be generated as a bagged, boosted and/or random forest ensemble.

The pre package is developed to provide useRs a completely R based implementation of the algorithm described by Friedman & Popescu (2008). However, note that pre is under development, and much work still needs to be done. See Fokkema, Smits, Kelderman, & Penninx (2015) for an application of the methods.

Examples
========

To get a first impression of how pre works, we will fit a prediction rule ensemble to predict Ozone levels using the airquality dataset:

``` r
library(pre)
complete <- complete.cases(airquality)
set.seed(42)
airq.ens <- pre(Ozone ~ ., data = airquality[complete, ])
```

We can print the resulting ensemble:

``` r
print(airq.ens)
#> 
#> Final ensemble with cv error within 1se of minimum: 
#>   lambda =  1.607146
#>   number of terms = 14
#>   mean cv error (se) = 299.798 (74.28342) 
#> 
#>          rule    coefficient                          description
#>   (Intercept)   62.658099734                                 <NA>
#>        rule72  -13.401881493              Wind > 5.7 & Temp <= 84
#>       rule216    8.166292702         Wind <= 10.3 & Solar.R > 148
#>       rule122    8.027236520                            Temp > 77
#>       rule213   -7.901556274              Wind > 5.1 & Temp <= 87
#>       rule201   -6.587690267  Wind > 5.7 & Temp <= 87 & Day <= 23
#>        rule25   -5.524545249              Wind > 6.3 & Temp <= 82
#>       rule179   -4.981266386              Wind > 5.7 & Temp <= 82
#>         rule3    4.927105788              Temp > 78 & Wind <= 6.3
#>       rule149   -3.427067580                Temp <= 87 & Wind > 8
#>       rule212    3.315627932                        Solar.R > 201
#>        rule89    2.588377937              Temp > 77 & Wind <= 8.6
#>        rule76   -2.112110981              Wind > 6.3 & Temp <= 84
#>       rule174   -1.315368661                           Wind > 6.9
#>       rule119   -0.009507402                Wind > 8 & Temp <= 76
```

We can plot a subsample of the rules (and or/linear terms) in the ensemble:

``` r
plot(airq.ens, max.terms.plot = 9, cex = .75)
```

![](inst/README-figures/README-unnamed-chunk-4-1.png)![](inst/README-figures/README-unnamed-chunk-4-2.png)

We can obtain the estimated coefficients of the ensemble:

``` r
head(coef(airq.ens))
#>            rule coefficient                         description
#> 194 (Intercept)   62.658100                                <NA>
#> 56       rule72  -13.401881             Wind > 5.7 & Temp <= 84
#> 170     rule216    8.166293        Wind <= 10.3 & Solar.R > 148
#> 95      rule122    8.027237                           Temp > 77
#> 167     rule213   -7.901556             Wind > 5.1 & Temp <= 87
#> 159     rule201   -6.587690 Wind > 5.7 & Temp <= 87 & Day <= 23
```

We can assess the importance of the predictor variables, and each of the rules and/or linear terms in the ensemble:

``` r
importance(airq.ens)
```

![](inst/README-figures/README-unnamed-chunk-6-1.png)

    #> $varimps
    #>   varname       imp
    #> 1    Temp 14.992106
    #> 2    Wind 13.550713
    #> 3 Solar.R  3.691361
    #> 4     Day  1.064621
    #> 
    #> $baseimps
    #>       rule                         description         imp   coefficient
    #> 1   rule72             Wind > 5.7 & Temp <= 84 6.098127160 -13.401881493
    #> 2  rule216        Wind <= 10.3 & Solar.R > 148 4.053274590   8.166292702
    #> 3  rule122                           Temp > 77 4.011974345   8.027236520
    #> 4  rule201 Wind > 5.7 & Temp <= 87 & Day <= 23 3.193863036  -6.587690267
    #> 5  rule213             Wind > 5.1 & Temp <= 87 3.108749351  -7.901556274
    #> 6   rule25             Wind > 6.3 & Temp <= 82 2.633249021  -5.524545249
    #> 7  rule179             Wind > 5.7 & Temp <= 82 2.342414876  -4.981266386
    #> 8    rule3             Temp > 78 & Wind <= 6.3 1.782470876   4.927105788
    #> 9  rule149               Temp <= 87 & Wind > 8 1.677079000  -3.427067580
    #> 10 rule212                       Solar.R > 201 1.664724095   3.315627932
    #> 11  rule89             Temp > 77 & Wind <= 8.6 1.198545921   2.588377937
    #> 12  rule76             Wind > 6.3 & Temp <= 84 0.985824758  -2.112110981
    #> 13 rule174                          Wind > 6.9 0.543944896  -1.315368661
    #> 14 rule119               Wind > 8 & Temp <= 76 0.004559346  -0.009507402
    #>           sd
    #> 1  0.4550202
    #> 2  0.4963421
    #> 3  0.4997952
    #> 4  0.4848229
    #> 5  0.3934351
    #> 6  0.4766454
    #> 7  0.4702449
    #> 8  0.3617683
    #> 9  0.4893627
    #> 10 0.5020841
    #> 11 0.4630490
    #> 12 0.4667486
    #> 13 0.4135304
    #> 14 0.4795575

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
#> [1] 365.9738
#> 
#> $MAE
#> [1] 13.77589
```

References
==========

Fokkema, M., Smits, N., Kelderman, H., & Penninx, B. W. (2015). Connecting clinical and actuarial prediction with rule-based methods. *Psychological Assessment*, *27*(2), 636.

Friedman, J. H., & Popescu, B. E. (2008). Predictive learning via rule ensembles. *The Annals of Applied Statistics*, *2*(3), 916–954. Retrieved from <http://www.jstor.org/stable/30245114>
