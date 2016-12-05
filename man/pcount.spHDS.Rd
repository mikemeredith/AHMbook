\name{pcount.spHDS}
\alias{pcount.spHDS}
\encoding{UTF-8}

\title{
Fit spatial hierarchical distance sampling model
}
\description{
Function fits spatial hierarchical distance sampling model; unmarked model fitting function introduced in Section 9.8.4. See \code{\link[unmarked]{pcount}} for details.
}
\usage{
pcount.spHDS(formula, data, K, mixture = c("P", "NB", "ZIP"), starts, method = "BFGS", se = TRUE, engine = c("C", "R"), ...)
}
\arguments{
  \item{formula}{
Double right-hand side formula describing covariates of detection and abundance, in that order
}
  \item{data}{
an unmarkedFramePCount object supplying data to the model.
}
  \item{K}{
Integer upper index of integration for N-mixture. This should be set high enough so that it does not affect the parameter estimates. Note that computation time will increase with K.
}
  \item{mixture}{
character specifying mixture: "P", "NB", or "ZIP".
}
  \item{starts}{
vector of starting values
}
  \item{method}{
Optimization method used by \code{\link{optim}}.
}
  \item{se}{
logical specifying whether or not to compute standard errors.
}
  \item{engine}{
Either "C" or "R" to use fast C++ code or native R code during the optimization.
}
  \item{\dots}{
Additional arguments to optim, such as lower and upper bounds
}
}

\value{
unmarkedFit object describing the model fit.
}
\references{
Kéry & Royle (2016) \emph{Applied Hierarachical Modeling in Ecology} Section 9.8.4
}
\author{
Kéry & Royle
}

\examples{
# TO DO
}
