% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FHme.R
\name{FHme}
\alias{FHme}
\title{Fay-Herriot Model with Measurement Error}
\usage{
FHme(
  formula,
  vardir,
  var.x,
  type.x = "witherror",
  MAXITER = 1000,
  PRECISION = 1e-04,
  data
)
}
\arguments{
\item{formula}{an object of class \code{\link[stats]{formula}} (or one that can be coerced to that class): a symbolic description of the model to be fitted. The variables included \code{formula} must have a length equal to the number of domains \code{m}. This formula can provide auxiliary variable either measured with error or without error or combination between them. If the auxiliary variable are combination between \code{noerror} and \code{witherror} variable, input all \code{witherror} variable first then \code{noerror} variable.}

\item{vardir}{vector containing the \code{m} sampling variances of direct estimators for each domain. The values must be sorted as the \code{Y}.}

\item{var.x}{vector containing mean squared error of \code{X} . The values must be sorted as the \code{X}. if you use optional \code{data}, input this parameter use \code{c("")}, example: \code{var.x = c("c1") or var.x = c("c1","c2")}.}

\item{type.x}{type of auxiliary variable used in the model. Either source measured with \code{noerror}, \code{witherror} and \code{mix}. Default value is \code{witherror}.}

\item{MAXITER}{maximum number of iterations allowed. Default value is \code{1000} iterations.}

\item{PRECISION}{convergence tolerance limit. Default value is \code{0.0001}.}

\item{data}{optional data frame containing the variables named in formula, vardir, and var.x.}
}
\value{
The function returns a list with the following objects:
\describe{
   \item{\code{eblup}}{vector with the values of the estimators for the domains.}
   \item{\code{fit}}{a list containing the following objects:}
   \itemize{
    \item \code{method} : type of fitting method.
    \item \code{convergence} : a logical value of convergence when calculating estimated beta and estimated random effects.
    \item \code{iterations} : number of iterations when calculating estimated beta and estimated random effects.
    \item \code{estcoef} : a data frame with the estimated model coefficient (\code{beta}) in the first column, their standard error (\code{std.error}) in the second column, the t-statistics (\code{t.statistics}) in the third column, and the p-values of the significance of each coefficient (\code{pvalue}) in the last column.
    \item \code{refvar} : a value of estimated random effects.
    \item \code{gamma} : vector with values of the estimated gamma for each domains.
    }
 }
}
\description{
This function gives the EBLUP estimator based on Fay-Herriot model with measurement error.
}
\details{
A formula has an implied intercept term. To remove this use either y ~ x - 1 or y ~ 0 + x. See \code{\link[stats]{formula}}  for more details of allowed formulae.
}
\examples{
data(dataME)
data(datamix)
sae.me <- FHme(formula = y ~ x.hat, vardir = vardir, var.x = c("var.x"), data = dataME)
sae.mix <- FHme(formula = y ~ x.hat1 + x.hat2 + x3 + x4,
            vardir = vardir, var.x = c("var.x1", "var.x2"), type.x = "mix", data = datamix)

}
\seealso{
\code{\link{mse_FHme}}
}
