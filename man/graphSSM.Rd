\name{graphSSM}
\alias{graphSSM}
\encoding{UTF-8}

% former name 'graph.ssm2'

\title{
Plot the results for a Gaussian state-space model
}
\description{
Plot the observed time-series and the estimated population trajectories for a multivariate Gaussian state-space model (SSM).
}
\usage{
graphSSM(ssm, C)
}
\arguments{
  \item{ssm}{
An object of class \code{jagsUI} with the results of fitting the SSM model.
}
  \item{C}{
The original count data: an individuals x year matrix of counts, no missing values.
}
}
\value{
No output, the function is called for plotting.
}
\references{
Kéry, M. & Royle, J.A. (2021) \emph{Applied Hierarchical Modeling in Ecology} Section 1.6

Kéry, M. & Schaub, M.  (2012) \emph{Bayesian population analysis using WinBUGS - a hierarchical perspective}, Academic Press.

Kéry, M. & Royle, J.A. (2021) \emph{Applied Hierarchical Modeling in Ecology} AHM2 - 1.6.
}
\author{
Marc Kéry & Andy Royle, adapted for multivariate SSMs from the SSM graphing function \code{graph.ssm} in Chapter 5 in Kéry & Schaub (2012).
}
