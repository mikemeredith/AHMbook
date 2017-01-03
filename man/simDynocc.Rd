\name{simDynocc}
\alias{simDynocc}
\encoding{UTF-8}
\title{
Simulate data under a non-spatial dynamic occupancy model
}
\description{
Function to simulate detection/nondetection data under a general dynamic site-occupancy model, including:

* annual variation in the probabilities of patch persistence, colonization and detection is specified by the bounds of a uniform distribution.

* one covariate is allowed to affect a parameter: a site covariate for psi1, site-by-year covariates for phi and gamma, and an observational covariate for p

* Additional detection heterogeneity at the site-, the survey, or the site-by-survey level, with the possibility of a temporal trend in this heterogeneity over the years. E.g., an annual trend in detection heterogeneity at the site or the survey level is specified by the first and second value, which correspond to the heterogeneity in the first and the last year Hence, range.sd.site = c(0, 1) will result in a linear trend in the magnitude of site heterogeneity in detection from 0 in the first year to 1 in the last year.

* Additional detection heterogeneity that varies over the season (= occasion) according to a quadratic effect of occasion number (to model phenology of an insect species for instance).

* Simulation of data under a BACI (before-after-control-impact) design, where some event happens in a given year and \emph{reduces} phi or gamma by a stated percentage (only reductions, no increases allowed !)

}
\usage{
simDynocc(nsite = 250, nrep = 3, nyear = 10, year.impact = 5,
  mean.psi1 = 0.4, beta.Xpsi1 = 0,
  range.phi = c(0.5, 1), impact.phi = 0, beta.Xphi = 0,
  range.gamma = c(0, 0.5), impact.gamma = 0, beta.Xgamma = 0,
  range.p = c(0.1, 0.9), beta.Xp = 0,
  range.beta1.season = c(0, 0), range.beta2.season = c(0, 0),
  range.sd.site = c(0, 0), range.sd.survey = c(0, 0),
  range.sd.site.survey = c(0, 0),
  show.plot = TRUE)
}
\arguments{
  \item{nsite}{
Number of sites
}
  \item{nrep}{
Number of replicate surveys within a year (= season)
}
  \item{nyear}{
Number of years (or 'seasons')
}
  \item{year.impact}{
Year when some impact happens (for BACI design)
}
  \item{mean.psi1}{
average occupancy probability in first year
}
  \item{beta.Xpsi1}{
coefficients of environmental covariate in probability of initial occupancy.
}
  \item{range.phi}{
bounds of uniform distribution from which annual persistence is drawn
}
  \item{impact.phi}{
effect in percent on annual phi (must be zero or negative, e.g., impact.phi = -20 for a 20\% reduction in phi)
}
  \item{beta.Xphi}{
coefficients of environmental covariate in probability of persistence.
}
  \item{range.gamma}{
bounds of uniform distribution from which annual colonisation is drawn
}
  \item{impact.gamma}{
effect in percent on annual gamma (must be zero or negative, e.g., impact.gamma = -20 for a 20\% reduction in gamma)
}  \item{beta.Xgamma}{
coefficients of environmental covariate in probability of colonization.
}
  \item{range.p}{
bounds of uniform distribution from which annual p is drawn
}
  \item{beta.Xp}{
coefficients of environmental covariate in probability of detection.
}
 \item{range.beta1.season}{
the range of the annual variation in the linear effect of season (i.e., of month 1-12) on the product of availability and detection linear and quadratic effect of season 
}
  \item{range.beta2.season}{
the same for the quadratic effect of season
}
  \item{range.sd.site}{
sd of normal distribution to model logit-normal noise in p at the site level in the first and the last year of the simulation; if the two values in the range are the same, a constant value is assumed over time, while if they are different, a linear trend is assumed over time.
}
  \item{range.sd.survey}{
sd of normal distribution to model logit-normal noise in p \emph{only} at the rep = 'survey' level, in the first and the last year; if they are different, a linear trend is assumed over time.
}
  \item{range.sd.site.survey}{
sd of normal distribution to model logit-normal noise in p at the site/year/rep = ‘survey’ level, in the first and the last year; if they are different, a linear trend is assumed over time.
}
  \item{show.plot}{
If TRUE, plots of results are displayed; set to FALSE if running many simulations.
}
}

\value{
A list with the values of the arguments input and the following additional elements:
  \item{impact }{a 0/1 vector of length nyear - 1 indicating if an impact applies to the interval}
  \item{phi.effect }{additive effect of impact on persistence, a vector of length nyear - 1}
  \item{gamma.effect }{additive effect of impact on colonisation, a vector of length nyear - 1}
  \item{beta1 }{linear effect of occasion on the product of availability and detection, a vector of length nyear}
  \item{beta2 }{quadratic effect of occasion on the product of availability and detection, a vector of length nyear}
  \item{mean.phi }{mean persistence for each interval before application of any BACI effect, a vector of length nyear - 1}
  \item{mean.gamma }{mean colonisation for each interval before application of any BACI effect, a vector of length nyear - 1}
  \item{mean.p }{mean detection probability for each year, a vector of length nyear}
  \item{psi }{annual occupancy for each site, a nsite x nyear matrix}
  \item{mean.psi }{average occupancy over sites, a vector of length nyear}
  \item{n.occ }{number of occupied sites, a vector of length nyear}
  \item{psi.fs }{finite-sample occupancy proportion, a vector of length nyear}
  \item{psi.app }{apparent occupancy over sites, a vector of length nyear}
  \item{z }{true occurrence state, a nsite x nyear matrix of 0/1}
  \item{phi }{persistence, a nsite x nyear-1 matrix}
  \item{gamma }{colonisation, a nsite x nyear-1 matrix}
  \item{p }{detection probability, a nsite x nrep x nyear array}
  \item{y }{the observed detection history, a nsite x nrep x nyear array}
  \item{Xpsi1 }{covariate affecting initial occupancy, a vector of length nsite}
  \item{Xphi }{covariate affecting persistence, a nsite x nyear matrix}
  \item{Xgamma }{covariate affecting colonisation, a nsite x nyear matrix}
  \item{Xp }{covariate affecting probability of detection, a nsite x nrep x nyear array}
  \item{eps3 }{the array of site/rep/year random effects, a nsite x nrep x nyear array}
}
\references{
Kéry & Royle (2017) \emph{Applied Hierarachical Modeling in Ecology} Chapter 16

}
\author{
Marc Kéry, 4 Dec 2014, modified 18-20 October 2016
}

\examples{
# Generate data with the default arguments and look at the structure:
tmp <- simDynocc()
str(tmp)
}
