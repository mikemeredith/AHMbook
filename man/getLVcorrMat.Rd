\name{getLVcorrMat}
\alias{getLVcorrMat}
\alias{getEcorrMat}
\encoding{UTF-8}

\title{
Compute results of a Joint Species Distribution Model (JSDM)
}
\description{
Compute pairwise residual correlation in occurrence or abundance in a Joint Species Distribution Model (JSDM) with latent variables.

\code{getLVcorrMat} computes the correlation matrix in residual occurrence or abundance in an latent-variable multi-species occupancy or N-mixture model which has been fit using MCMC (eg, in JAGS or WinBUGS). See Tobler et al. (2019) and Kéry & Royle (2021) chapter 8.

\code{getEcorrMat} computes the correlation matrix for the coefficients of the observed environmental covariates.
}
\usage{
getLVcorrMat(lv.coef, type=c("occupancy", "Nmix"), stat=mean)

getEcorrMat(beta, stat=mean)
}
\arguments{
  \item{lv.coef}{
MCMC chains for the coefficients of the latent variables, typically from the \code{sims.list} of a model fit; an iterations x species x latent variables array.
}
  \item{beta}{
MCMC chains for the coefficients of the environment variables (excluding the intercept), typically from the \code{sims.list} of a model fit; an iterations x species x environmental variables array.
}
  \item{type}{
Indication of whether the model fit was for occupancy data (with a probit link) or an N-mixture model based on count data.
}
  \item{stat}{
The function used to summarize the MCMC chains for the correlations; if \code{stat = NULL}. the full array with the MCMC chains is returned
}
}
\value{
The relevant correlation matrix, species x species; if \code{stat = NULL}, an array, iterations x species x species, with the MCMC chains.
}
\references{
Tobler, M. et al. (2019) Joint species distribution models with species correlations and imperfect detection. \emph{Ecology}, 100(8), e02754.

Kéry, M. & Royle, J.A. (2021) \emph{Applied Hierarchical Modeling in Ecology} AHM2 - 8.4.2 and 8.5.4.
}
\author{
Mathias Tobler.
}

