\name{spline.prep}
\alias{spline.prep}
\encoding{UTF-8}
\title{
Prepare input for BUGS model when fitting a spline for a covariate
}
\description{
Function chooses knots and creates design matrices for fixed and random-effects parts of a spline model for a chosen covariate. Based on code by Crainiceanu et al. (2005) and Zuur et al. (2012). Allows you to choose number of knots or else uses it by the rule given in Crainiceanu et al. (2005). Prepares fixed part of covariate as a quadratic polynomial.
}
\usage{
spline.prep(cov, nknot = NA)
}
\arguments{
  \item{cov}{
the covariate, a numeric vector
}
  \item{nknot}{
optional, number of knots
}
}
\value{
A list with the following elements:
  \item{cov }{The input covariate}
  \item{knots }{The values of the knots}
  \item{X }{The fixed-effects design matrix}
  \item{Z }{The random-effects design matrix}
}
\references{
Crainiceanu, C.M., Ruppert, D., & Wand, M.P. (2005) Bayesian analysis for penalized spline regression using WinBUGS. \emph{Journal of Statistical Software}, 14.

Zuur, A.F., Saveliev, A.A., Ieno, E.N. (2012) \emph{Zero-inflated Models and Generalized Linear Mixed Models with R}. Highlands Statistics

Kéry, M. & Royle, J.A. (2016) \emph{Applied Hierarchical Modeling in Ecology} AHM1 - 10.14.
}
\author{
Marc Kéry & Andy Royle
}

