\name{simDCM}
\alias{simDCM}
\encoding{UTF-8}

\title{
Simulate detection/nondetection data under a general dynamic community (site-occupancy) model
}
\description{
Function to simulate detection/nondetection data under a general dynamic community (= dynamic, multi-species site-occupancy) model, including:

* annual variation in the probabilities of patch persistence, colonization and detection is specified by the bounds of a uniform distribution.

* species heterogeneity around the means is specified by the SD of a normal distribution and expressed on the logit scale

* one covariate is allowed to a parameter (site covariate for psi1, site-year covariate for phi and gamma and site-year-rep covariate for p). Each covariate is allowed to differ among species again according to a logit-normal model of heterogeneity.

* additional detection heterogeneity at the site- or the occasion level, with the possibility of a temporal trend in this heterogeneity over years. E.g., an annual trend in detection heterogeneity at the site or the occasion level is specified by the value in the first and the last year. Hence, range.sd.site = c(0, 1) will result in a linear trend in the magnitude of site-level heterogeneity in detection from 0 in the first year to 1 in the last year, with interpolation for the years in between.

* additional detection heterogeneity that among occasions according to a quadratic effect of occasion number (to model phenology of an insect species for instance).

These last two types of detection heterogeneity are not (yet) allowed to be species-specific.
}
\usage{
simDCM(nspec = 50, nsite = 100, nrep = 3, nyear = 10, 
  mean.psi1 = 0.4, sig.lpsi1 = 1, mu.beta.lpsi1 = 0, sig.beta.lpsi1 = 0,
  range.mean.phi = c(0.8, 0.8), sig.lphi = 1, mu.beta.lphi = 0,
  sig.beta.lphi = 0, range.mean.gamma = c(0.2, 0.2), sig.lgamma = 1,
  mu.beta.lgamma = 0, sig.beta.lgamma = 0, range.mean.p = c(0.5, 0.5),
  sig.lp = 1, mu.beta.lp = 0, sig.beta.lp = 0, range.beta1.season = c(0, 0),
  range.beta2.season = c(0, 0), range.sd.site = c(0, 0),
  range.sd.survey = c(0, 0), show.plot = TRUE)
}

\arguments{
  \item{nspec}{
number of species (typically called N in AHM book)
}
  \item{nsite}{
number of sites (M)
}
  \item{nrep}{
number of replicate occasions within a year (J)
}
  \item{nyear}{
number of years (T)
}
  \item{mean.psi1}{
average (across all species in the community) of the intercept of occupancy probability in first year
}
  \item{sig.lpsi1}{
sd of the normal distribution from which species-specific occupancy intercepts are drawn (centered on logit(mean.psi1)), on logit scale
}
  \item{mu.beta.lpsi1}{
community mean of the coefficients of the covariate in probabilities of initial occupancy: the probability-scale mean of the normal distribution from which the species-specific coefficients are drawn.
}
  \item{sig.beta.lpsi1}{
sd of the normal distribution from which species-specific slopes are drawn (centered on mu.beta.lpsi1)
}
  \item{range.mean.phi}{
bounds of uniform distribution from which the average (across species) annual intercept of persistence is drawn
}
  \item{sig.lphi}{
sd of the normal distribution from which species-specific persistence intercepts are drawn (centered on logit(mean.phi), which are year-specific), on logit scale
}
  \item{mu.beta.lphi}{
community mean of the coefficients of the covariate in probabilities of persistence: the probability-scale mean of the normal distribution from which the species-specific coefficients are drawn.
}
  \item{sig.beta.lphi}{
sd of the normal distribution from which species-specific persistence slopes are drawn (centered on mu.beta.lphi)
}
  \item{range.mean.gamma}{
bounds of uniform distribution from which the average (across species) annual intercept of colonisation is drawn
}
  \item{sig.lgamma}{
sd of the normal distribution from which species-specific colonization intercepts are drawn (centered on logit(mean.gamma), which are year-specific), on logit scale
}
  \item{mu.beta.lgamma}{
community mean of the coefficients of the covariate in probabilities of colonization: the probability-scale mean of the normal distibution from which the species-specific coefficients are drawn.
}
  \item{sig.beta.lgamma}{
sd of the normal distribution from which species-specific colonization slopes are drawn (centered on mu.beta.lgamma)
}
  \item{range.mean.p}{
bounds of uniform distribution from which the average (across species) annual intercept of p is drawn
}
  \item{sig.lp}{
sd of the normal distribution from which species-specific detection intercepts are drawn (centered on logit(mean.p), which are year-specific), on logit scale
}
  \item{mu.beta.lp}{
community mean of the coefficients of the covariate in probabilities of detection: the probability-scale mean of the normal distibution from which the species-specific coefficients are drawn.
}
  \item{sig.beta.lp}{
sd of the normal distribution from which species-specific detection slopes are drawn (centered on mu.beta.lp)
}
  \item{range.beta1.season}{
the range of the annual variation in the linear effect of season (i.e., of month 1-12 if nrep = 12) on the product of availability and detection 
}
  \item{range.beta2.season}{
the range of the annual variation in the quadratic effect of season (i.e., of month 1-12 if nrep = 12) on the product of availability and detection 
}
  \item{range.sd.site}{
sd of normal distribution to model logit-normal noise in p at the site level in the first and the last year of the simulation, with values for internediate years interpolated linearly
}
  \item{range.sd.survey}{
sd of normal distribution to model logit-normal noise in p at the occasion level, in the first and the last year, with values for internediate years interpolated linearly
}
  \item{show.plot}{
if TRUE, plots are produced. Set this to FALSE when running many simulations.
}
}

\value{
A list with the values of the input arguments and the following additional elements:
 \item{Xpsi1 }{Site covariate for psi1, a nsite x 1 matrix}
 \item{Xphi }{Yearly-site covariate for phi, a nsite x nyear matrix}
 \item{Xgamma }{Yearly-site covariate for gamma, a nsite x nyear matrix}
 \item{Xp }{Observation covariate for p, a nsite x nrep x nyear array}
 \item{beta0.lpsi }{initial (logit-scale) occupancy intercept for each species in the community, a vector of length nspec}
 \item{beta1.lpsi }{initial (log-scale) occupancy slope on Xpsi1 for each species in the community, a vector of length nspec}
 \item{psi }{occupancy probability per site, year and species, a nsite x nyear x nspec array}
 \item{mean.phi }{mean persistence (across species) intercept for each interval, a vector of length (nyear - 1)}
 \item{mean.gamma }{mean colonisation (across species) intercept for each interval, a vector of length (nyear - 1)}
 \item{eps.lphi }{additive species effects in logit(phi) intercept, a vector of length nspec}
 \item{eps.lgamma }{additive species effects in logit(gamma) intercept, a vector of length nspec}
 \item{beta0.lphi }{logit-scale persistence intercepts for each species in community, a nspec x (nyear - 1) matrix}
 \item{beta0.lgamma }{logit scale colonisation intercepts for each species in the community, a nspec x (nyear - 1) matrix}
 \item{beta1.lphi }{slope of logit(phi) on Xphi for each species in the community, a vector of length nspec}
 \item{beta1.lgamma }{slope of logit(gamma) on Xgamma for each species in the community, a vector of length nspec}
 \item{phi }{probability of persistence for each site, yearly interval and species, a nsite x (nyear-1) x nspec array}
 \item{gamma }{probability of colonisation for each site, yearly interval and species, a nsite x (nyear-1) x nspec array}
 \item{mean.p }{mean detection (across species) intercept for each year, a vector of length nyear}
 \item{eps.lp }{additive species effects in logit(p) intercept, a vector of length nspec}
 \item{beta0.lp }{species- and site-specific intercepts in the linear predictor for p, a nspec x nyear matrix}
 \item{beta1.lp }{species specific slopes of logit(p) on Xp, a vector of length nspec}
 \item{beta1 }{linear effect of the occasion number on detection probability, a vector of length nyear}
 \item{beta2 }{quadratic effect of the occasion number on detection probability, a vector of length nyear}
 \item{sd.site }{standard deviation of the zero-mean normal distribution from which additional, site-level detection heterogeneity is simulated, a vector of length nyear}
 \item{sd.survey }{standard deviation of the zero-mean normal distribution from which additional, occasion-level detection heterogeneity is simulated,, a vector of length nyear}
 \item{eps1 }{additive site random effects tht generate unstructured site-level detection heterogeneity, a vector of length nsite}
 \item{eps2 }{additive occasion random effects tht generate unstructured site-level detection heterogeneity, a vector of length nrep}
 \item{n.occ }{Number of occupied sites, a nyear x nspec matrix}
 \item{psi.fs }{Finite-sample occupancy proportion, a nyear x nspec matrix}
 \item{mean.psi }{Average psi over sites, a nyear x nspec matrix}
 \item{z.obs }{Observed value of z matrix, a nsite x nyear x nspec array}
 \item{n.occ.obs }{Observed number of occupied sites, a nyear x nspec matrix}
 \item{psi.obs }{Observed occupancy (finite sample), a nyear x nspec matrix}
 \item{nyear.pres }{Number of years when species present, a vector of length nspec}
 \item{nspec.pres }{Number of species ever present, scalar}
 \item{nyear.det }{Number of years when species detected, a vector of length nspec}
 \item{nspec.det }{Number of species ever detected, scalar}
 \item{z }{True value of z matrix (ie, presence/absence), a nsite x nyear x nspec array}
 \item{p }{Probability of detection, a nsite x nrep x nyear x nspec array}
 \item{y }{Observed detection history, a nsite x nrep x nyear x nspec array of 0/1}
}
\references{
Kéry & Royle (2017) \emph{Applied Hierarachical Modeling in Ecology} Section 17.2
}
\author{
Marc Kéry, 28 Nov 2016
}

\examples{
# Simulate a data set with the default arguments and look at the structure of the output:
tmp <- simDCM()
str(tmp)
}
