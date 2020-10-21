\name{simPOP}
\alias{simPOP}
\encoding{UTF-8}

\title{
Simulate data for a demographic state-space model
}
\description{
Simulate multiple time-series of counts under a pure Markov model (with exponential population model) or under an extended Markov model with exponential-plus-random-immigration population model; see Sollmann et al.(2015). Default is Markov model, setting \code{sd.rho} to a value greater than 0 changes to extended Markov and sets the amount of random immigration.
}
\usage{
simPOP(M = 100, T = 10, mean.lam = 3, beta.lam = 0, sd.log.lam = 0,
  mean.gamma = 1.0, beta.gamma = 0, sd.log.gamma.site = 0,
  sd.log.gamma.time = 0, sd.log.gamma.survey = 0, sd.rho = 0,
  mean.p = 0.6, beta.p = 0, sd.logit.p.site = 0, sd.logit.p.time = 0,
  sd.logit.p.survey = 0, show.plot = TRUE)
}
\arguments{
  \item{M}{
The number of sites.
}
  \item{T}{
The number of years.
}
  \item{mean.lam}{
The mean abundance for year 1.
}
  \item{beta.lam}{
The covariate coefficient for lambda.
}
  \item{sd.log.lam}{
The over-dispersion in lambda.
}
  \item{mean.gamma}{
The mean population growth rate.
}
  \item{beta.gamma}{
The covariate coefficient for gamma.
}
  \item{sd.log.gamma.site}{
SD of random site effects for gamma.
}
  \item{sd.log.gamma.time}{
SD of random time effects for gamma.
}
  \item{sd.log.gamma.survey}{
SD of random survey (site+time) effects for gamma.
}
  \item{sd.rho}{
The random immigration term.
}
  \item{mean.p}{
The mean detection probability.
}
  \item{beta.p}{
The covariate coefficient for p.
}
  \item{sd.logit.p.site}{
SD of random site effects for p on the logit scale.
}
  \item{sd.logit.p.time}{
SD of random time effects for p on the logit scale.
}
  \item{sd.logit.p.survey}{
SD of random survey (site+time) effects for p on the logit scale.
}
  \item{show.plot}{
Choose whether to show plots or not. Set to FALSE when using function in simulations.
}
}
\value{
A list with the values of the arguments entered and the following additional elements:

\item{Xsite1 }{M vector, site covariate affecting initial abundance (lambda).}
\item{Xsiteyear1 }{M x T matrix, yearly site covariate affecting recruitment (gamma).}
\item{Xsiteyear2 }{M x T matrix, yearly site covariate affecting detection (p).}
\item{eps.N }{M vector, site over-dispersion at t = 1.}
\item{lambda }{M vector, abundance in year 1.}
\item{eps.gamma.site }{M vector, random site effect for gamma.}
\item{eps.gamma.time }{T vector, random time effect for gamma.}
\item{eps.gamma.survey }{M x T matrix, random survey effect for gamma.}
\item{gamma }{M x T matrix, population growth rate.}
\item{rho }{(T-1) vector, immigration rate.}
\item{eps.p.site }{M vector, random site effect for detection.}
\item{eps.p.time }{T vector, random time effect for detection.}
\item{eps.p.survey }{M x T matrix, random survey effect for detection.}
\item{p }{M x T matrix, detection probability.}
\item{N }{M x T matrix, true population.}
\item{C }{M x T matrix, simulated counts.}
\item{zeroNyears }{scalar, sum(N == 0).}
\item{Nextinct }{scalar, number of sites where N == 0 at time T.}
\item{extrate }{scalar, proportion of sites where N == 0 at time T.}
\item{sumN }{T vector, total population in each year.}
\item{gammaX }{(T-1) vector, realized population growth rate.}

}
\references{
Sollmann, R. et al. (2015) An open-population hierarchical distance sampling model. \emph{Ecology}, 96, 325-331.

Kéry, M. & Royle, J.A. (2021) \emph{Applied Hierarchical Modeling in Ecology} AHM2 - 1.7.1.
}
\author{
Marc Kéry & Andy Royle
}

\examples{
# Run with the default arguments and look at the structure of the output
set.seed(123)
tmp <- simPOP()
str(tmp)
head(tmp$C)


}
