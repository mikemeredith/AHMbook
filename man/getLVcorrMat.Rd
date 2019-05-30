\name{getLVcorrMat}
\alias{getLVcorrMat}
\encoding{UTF-8}

\title{
Compute results of a species correlation model
}
\description{
Computes the correlation matrix in residual occurrence in an latent-variable multi-species occupancy or N-mixture model which has been fit using MCMC (eg, in JAGS or WinBUGS). See Tobler et al. (2019) and Kéry & Royle (2019) chapter 19.
}
\usage{
getLVcorrMat(lv.coef, type=c("occupancy", "Nmix"), stat=mean)
}
\arguments{
  \item{lv.coef}{
MCMC chains for the coefficients of the latent variables, typically from the \code{sims.list} of a model fit; an iterations x species x latent variables array.
}
  \item{type}{
Indication of whether the model fit was for occupancy data (with a probit link) or an N-mixture model based on count data.
}
  \item{stat}{
The function used to summarise the MCMC chains for the correlations.
}
}
\value{
The residual correlation matrix, species x species.
}
\references{
Tobler et al. (2019) Joint species distribution models wtih species correlations and imperfect detection. \emph{Ecology} (accepted) doi: 10.1002/ecy.2754

Kéry & Royle (2019) \emph{Applied Hierarchical Modeling in Ecology} Section 19.5.1
}
\author{
Mathias Tobler.
}

