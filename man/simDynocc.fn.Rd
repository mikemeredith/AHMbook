\name{simDynocc.fn}
\alias{simDynocc.fn}
\encoding{UTF-8}
\title{
Simulate data under a non-spatial dynamic occupancy model
}
\description{
Function to simulate detection/nondetection data under a very general dynamic site-occ model, including:

* annual variation in the probabilities of patch persistence, colonization and detection is specified by the bounds of a uniform distribution.

* one site-specific (???) covariate each is allowed to affect every parameter

* autologistic effects in persistence and colonization probability can be chosen, which fits a logistic regression of these parameters on the proportion of occupied neighbouring cells (in a queen's or 2nd order neighbourhood)

* Additional detection heterogeneity at the site- or the survey level, with the possibility of a temporal trend in this heterogeneity. E.g., an annual trend in detection heterogeneity at the site or the survey level is specified by the value in the first and the last year. Hence, range.sd.site = c(0, 1) will result in a linear trend in the magnitude of site heterogeneity in detection from 0 in the first year to 1 in the last year.
}
\usage{
simDynocc.fn(M = 250, J = 3, K = 10, mean.psi1 = 0.4, beta.Xpsi1 = 0,
 range.phi = c(0.5, 1), beta.Xphi = 0, range.gamma = c(0, 0.5), 
 beta.Xgamma = 0, range.p = c(0.1, 0.9), beta.Xp = 0, range.sd.site = c(0, 0),
 range.sd.survey = c(0, 0), seed = Sys.time(), show.plot = TRUE)
}
\arguments{
  \item{M}{
Number of sites
}
  \item{J}{
Number of replicate surveys within a year (= season)
}
  \item{K}{
Number of years (or 'seasons')
}
  \item{mean.psi1}{
average occupancy probability in first year
}
  \item{beta.Xpsi1}{
coefficients of environmental covariates in probability of initial occupancy.
}
  \item{range.phi}{
bounds of uniform distribution from which annual persistence is drawn
}
  \item{beta.Xphi}{
coefficients of environmental covariates in probability of persistence.
}
  \item{range.gamma}{
bounds of uniform distribution from which annual colonisation is drawn
}
  \item{beta.Xgamma}{
coefficients of environmental covariates in probability of colonization.
}
  \item{range.p}{
bounds of uniform distribution from which annual p is drawn
}
  \item{beta.Xp}{
coefficients of environmental covariates in probability of detection.
}
  \item{range.sd.site}{
sd of normal distribution to model logit-normal noise in p at the site level in the first and the last year of the simulation
}
  \item{range.sd.survey}{
sd of normal distribution to model logit-normal noise in p at the site/year/rep = 'survey' level, in the first and the last year
}
  \item{seed}{
a seed for random number generation
}
  \item{show.plot}{
If TRUE, plots of results are displayed; set to FALSE if running many simulations.
}
}

\value{
A list with the values of the arguments input and the following additional elements:
  \item{mean.phi }{mean persistence for each interval, a vector of length K - 1}
  \item{mean.gamma }{mean colonisation for each interval, a vector of length K - 1}
  \item{mean.p }{mean detection probability for each year, a vector of length K}
  \item{psi }{annual occupancy for each site, a M x K matrix}
  \item{mean.psi }{average occupancy over sites, a vector of length K}
  \item{psi.app }{apparent occupancy over sites, a vector of length K}
  \item{z }{true occurence state, a M x K matrix of 0/1}
  \item{phi }{true persistence, a M x K-1 matrix}
  \item{gamma }{true colonisation, a M x K-1 matrix}
  \item{p }{true detection probability, a M x J x K array}
  \item{y }{detection history, a M x J x K array}
}
\references{
Kéry & Royle (2016) \emph{Applied Hierarachical Modeling in Ecology} Section 15

}
\author{
Marc Kéry, 4 Dec 2014
}

\examples{
# Generate data with the default arguments and look at the structure:
tmp <- simDynocc.fn()
str(tmp)
}
