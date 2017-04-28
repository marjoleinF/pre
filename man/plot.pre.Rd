% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pre.R
\name{plot.pre}
\alias{plot.pre}
\title{Plot method for class pre}
\usage{
\method{plot}{pre}(x, penalty.par.val = "lambda.1se", linear.terms = TRUE,
  nterms = NULL, max.terms.plot = 16, ask = FALSE, exit.label = "0",
  ...)
}
\arguments{
\item{x}{an object of class \code{\link{pre}}.}

\item{penalty.par.val}{character. Which value of the penalty parameter
criterion should be used? The value yielding minimum cv error
(\code{"lambda.min"}) or penalty parameter yielding error within 1 standard
error of minimum cv error ("\code{lambda.1se}")? Alternatively, a numeric 
value may be specified, corresponding to one of the values of lambda in the 
sequence used by glmnet, for which estimated cv error can be inspected by 
running \code{x$glmnet.fit} and \code{plot(x$glmnet.fit)}.}

\item{linear.terms}{logical. Should linear terms be included in the plot?}

\item{nterms}{numeric. The total number of terms (or rules, if 
\code{linear.terms = FALSE}) to be plotted. Default is \code{NULL}, 
resulting in all terms of the final ensemble to be plotted.}

\item{max.terms.plot}{numeric. The maximum number of terms per plot. Rules 
are plotted in a square pattern, so \code{is.integer(sqrt(max.terms.plot))} 
should return \code{TRUE}, otherwise max.terms.plot will be set to the next 
higher value which returns true. The default \code{max.terms.plot = 16} 
results in max. 4x4 rules per plot. If the number of terms exceeds the value 
specified for max.rules.plot, multiple pages of plots will be created.}

\item{ask}{logical. Should user be prompted before starting a new page of
plots?}

\item{exit.label}{character string. What label should be printed in nodes to 
which the rule does not apply (``exit nodes'')?}

\item{...}{Currently not used.}
}
\description{
\code{plot.pre} creates one or more plots depicting the rules in the final
ensemble as simple decision trees.
}
\examples{
\donttest{
 set.seed(42)
 airq.ens <- pre(Ozone ~ ., data=airquality[complete.cases(airquality),])
 plot(airq.ens)}
}
