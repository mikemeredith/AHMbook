\name{simDCM}
\alias{simDCM}
\encoding{UTF-8}

\title{
Simulate detection/nondetection data under a general dynamic community (site-occupancy) model
}
\description{
Function to simulate detection/nondetection data under a general dynamic community (= dynamic, multi-species site-occupancy) model, including:

* annual variation in the probabilities of patch persistence, colonization and detection is specified by the bounds of a uniform distribution.

* species heterogeneity around the community means is specified by the SD of a normal distribution and expressed on the logit scale

* one covariate is allowed per parameter (site covariate for psi1, site-year covariate for phi and gamma and site-year-rep covariate for p). Each covariate is allowed to differ among species again according to a logit-normal model of heterogeneity.

* additional detection heterogeneity at the site- or the occasion level, with the possibility of a temporal trend in this heterogeneity over years. E.g., an annual trend in detection heterogeneity at the site or the occasion level is specified by the value in the first and the last year. Hence, trend.sd.site = c(0, 1) will result in a linear trend in the magnitude of site-level heterogeneity in detection from 0 in the first year to 1 in the last year, with interpolation for the years in between.

* additional detection heterogeneity that among occasions according to a quadratic effect of occasion number (to model the typical 'phenological curve' of an insect species for instance).

These last two types of detection heterogeneity are not (yet) allowed to be species-specific.
}
\usage{
simDCM(nspecies = 50, nsites = 100, nsurveys = 3, nyears = 10,
  mean.psi1 = 0.4, sig.lpsi1 = 1, mu.beta.lpsi1 = 0, sig.beta.lpsi1 = 0,
  range.mean.phi = c(0.8, 0.8), sig.lphi = 1, mu.beta.lphi = 0,
  sig.beta.lphi = 0, range.mean.gamma = c(0.2, 0.2), sig.lgamma = 1,
  mu.beta.lgamma = 0, sig.beta.lgamma = 0, range.mean.p = c(0.5, 0.5),
  sig.lp = 1, mu.beta.lp = 0, sig.beta.lp = 0, range.beta1.survey = c(0, 0),
  range.beta2.survey = c(0, 0), trend.sd.site = c(0, 0),
  trend.sd.survey = c(0, 0), show.plot = TRUE, verbose = TRUE)
}

\arguments{
  \item{nspecies}{
number of species (typically called N in AHM book)
}
  \item{nsites}{
number of sites (M).
}
  \item{nsurveys}{
number of replicate surveys within a year (J).
}
  \item{nyears}{
number of years (T).
}
  \item{mean.psi1}{
average (across all species in the community) of the intercept of occupancy probability in first year.
}
  \item{sig.lpsi1}{
sd of the normal distribution from which species-specific occupancy intercepts are drawn (centered on logit(mean.psi1)), on logit scale.
}
  \item{mu.beta.lpsi1}{
community mean of the coefficients of the covariate in probabilities of initial occupancy: the probability-scale mean of the normal distribution from which the species-specific coefficients are drawn.
}
  \item{sig.beta.lpsi1}{
sd of the normal distribution from which species-specific slopes are drawn (centered on mu.beta.lpsi1).
}
  \item{range.mean.phi}{
bounds of uniform distribution from which the average (across species) annual intercept of persistence is drawn.
}
  \item{sig.lphi}{
sd of the normal distribution from which species-specific persistence intercepts are drawn (centered on logit(mean.phi), which are year-specific), on logit scale.
}
  \item{mu.beta.lphi}{
community mean of the coefficients of the covariate in probabilities of persistence: the probability-scale mean of the normal distribution from which the species-specific coefficients are drawn.
}
  \item{sig.beta.lphi}{
sd of the normal distribution from which species-specific persistence slopes are drawn (centered on mu.beta.lphi).
}
  \item{range.mean.gamma}{
bounds of uniform distribution from which the average (across species) annual intercept of colonization is drawn.
}
  \item{sig.lgamma}{
sd of the normal distribution from which species-specific colonization intercepts are drawn (centered on logit(mean.gamma), which are year-specific), on logit scale.
}
  \item{mu.beta.lgamma}{
community mean of the coefficients of the covariate in probabilities of colonization: the probability-scale mean of the normal distribution from which the species-specific coefficients are drawn.
}
  \item{sig.beta.lgamma}{
sd of the normal distribution from which species-specific colonization slopes are drawn (centered on mu.beta.lgamma).
}
  \item{range.mean.p}{
bounds of uniform distribution from which the average (across species) annual intercept of p is drawn.
}
  \item{sig.lp}{
sd of the normal distribution from which species-specific detection intercepts are drawn (centered on logit(mean.p), which are year-specific), on logit scale.
}
  \item{mu.beta.lp}{
community mean of the coefficients of the covariate in probabilities of detection: the probability-scale mean of the normal distribution from which the species-specific coefficients are drawn.
}
  \item{sig.beta.lp}{
sd of the normal distribution from which species-specific detection slopes are drawn (centered on mu.beta.lp).
}
  \item{range.beta1.survey}{
the range of the annual variation in the linear effect of survey (i.e., of surveys 1-12 if nsurveys = 12) on the product of availability and detection.
}
  \item{range.beta2.survey}{
the range of the annual variation in the quadratic effect of survey (i.e., of surveys 1-12 if nsurveys = 12) on the product of availability and detection.
}
  \item{trend.sd.site}{
sd of normal distribution to model logit-normal noise in p at the site level in the first and the last year of the simulation, with values for intermediate years interpolated linearly.
}
  \item{trend.sd.survey}{
sd of normal distribution to model logit-normal noise in p at the occasion level, in the first and the last year, with values for intermediate years interpolated linearly.
}
  \item{show.plot}{
if TRUE, plots are produced. Set this to FALSE when running simulations.
}
  \item{verbose}{
if TRUE, output will be written to the console.
}
}

\value{
A list with the values of the input arguments and the following additional elements:
 \item{Xpsi1 }{Site covariate for psi1, a nsites x 1 matrix}
 \item{Xphi }{Yearly-site covariate for phi, a nsites x nyears matrix}
 \item{Xgamma }{Yearly-site covariate for gamma, a nsites x nyears matrix}
 \item{Xp }{Observation covariate for p, a nsites x nsurveys x nyears array}
 \item{beta0.lpsi }{initial (logit-scale) occupancy intercept for each species in the community, a vector of length nspecies}
 \item{beta1.lpsi }{initial (log-scale) occupancy slope on Xpsi1 for each species in the community, a vector of length nspecies}
 \item{psi }{occupancy probability per site, year and species, a nsites x nyears x nspecies array}
 \item{mean.phi }{mean persistence (across species) intercept for each interval, a vector of length (nyears - 1)}
 \item{mean.gamma }{mean colonization (across species) intercept for each interval, a vector of length (nyears - 1)}
 \item{eps.lphi }{additive species effects in logit(phi) intercept, a vector of length nspecies}
 \item{eps.lgamma }{additive species effects in logit(gamma) intercept, a vector of length nspecies}
 \item{beta0.lphi }{logit-scale persistence intercepts for each species in community, a nspecies x (nyears - 1) matrix}
 \item{beta0.lgamma }{logit scale colonization intercepts for each species in the community, a nspecies x (nyears - 1) matrix}
 \item{beta1.lphi }{slope of logit(phi) on Xphi for each species in the community, a vector of length nspecies}
 \item{beta1.lgamma }{slope of logit(gamma) on Xgamma for each species in the community, a vector of length nspecies}
 \item{phi }{probability of persistence for each site, yearly interval and species, a nsites x (nyears-1) x nspecies array}
 \item{gamma }{probability of colonization for each site, yearly interval and species, a nsites x (nyears-1) x nspecies array}
 \item{mean.p }{mean detection (across species) intercept for each year, a vector of length nyears}
 \item{eps.lp }{additive species effects in logit(p) intercept, a vector of length nspecies}
 \item{beta0.lp }{species- and site-specific intercepts in the linear predictor for p, a nspecies x nyears matrix}
 \item{beta1.lp }{species specific slopes of logit(p) on Xp, a vector of length nspecies}
 \item{beta1 }{linear effect of the occasion number on detection probability, a vector of length nyears}
 \item{beta2 }{quadratic effect of the occasion number on detection probability, a vector of length nyears}
 \item{sd.site }{standard deviation of the zero-mean normal distribution from which additional, site-level detection heterogeneity is simulated, a vector of length nyears}
 \item{sd.survey }{standard deviation of the zero-mean normal distribution from which additional, occasion-level detection heterogeneity is simulated,, a vector of length nyears}
 \item{eps1 }{additive site random effects tht generate unstructured site-level detection heterogeneity, a vector of length nsites}
 \item{eps2 }{additive occasion random effects tht generate unstructured site-level detection heterogeneity, a vector of length nsurveys}
 \item{n.occ }{Number of occupied sites, a nyears x nspecies matrix}
 \item{psi.fs }{Finite-sample occupancy proportion, a nyears x nspecies matrix}
 \item{mean.psi }{Average psi over sites, a nyears x nspecies matrix}
 \item{z.obs }{Observed value of z matrix, a nsites x nyears x nspecies array}
 \item{n.occ.obs }{Observed number of occupied sites, a nyears x nspecies matrix}
 \item{psi.obs }{Observed occupancy (finite sample), a nyears x nspecies matrix}
 \item{nyears.pres }{Number of years when species present, a vector of length nspecies}
 \item{nspecies.pres }{Number of species ever present, scalar}
 \item{nyears.det }{Number of years when species detected, a vector of length nspecies}
 \item{nspecies.det }{Number of species ever detected, scalar}
 \item{z }{True value of z matrix (ie, presence/absence), a nsites x nyears x nspecies array}
 \item{p }{Probability of detection, a nsites x nsurveys x nyears x nspecies array}
 \item{y }{Observed detection history, a nsites x nsurveys x nyears x nspecies array of 0/1}
}
\references{
Kéry, M. & Royle, J.A. (2021) \emph{Applied Hierarchical Modeling in Ecology} AHM2 - 5.2.
}
\author{
Marc Kéry & Andy Royle
}

\examples{
# Simulate a data set with the default arguments and look at the structure of the output:
tmp <- simDCM()
str(tmp)

\donttest{
# Default arguments, without plots
str(data <- simDCM(show.plot = FALSE))
# More examples:
str(data <- simDCM(nspecies = 200)) # More species (looks great)
str(data <- simDCM(nspecies = 1))   # A single species (ha, works !)
str(data <- simDCM(nsites = 267))   # More sites
str(data <- simDCM(nsites = 1))     # A single site
str(data <- simDCM(nsurveys = 10))  # More visits
str(data <- simDCM(nyears = 25))    # More years
str(data <- simDCM(nyears = 2))     # Just two years
try(data <- simDCM(nyears = 1))     # A single year ... error

# No species heterogeneity in parameters of initial occupancy
str(data <- simDCM(sig.lpsi1 = 0, sig.beta.lpsi1 = 0))
# No species heterogeneity in parameters of persistence
str(data <- simDCM(sig.lphi = 0, sig.beta.lphi = 0))
# No species heterogeneity in parameters of colonisation
str(data <- simDCM(sig.lgamma = 0, sig.beta.lgamma = 0))
# No species heterogeneity in parameters of detection
str(data <- simDCM(sig.lp = 0, sig.beta.lp = 0))
# No annual variation in rates
str(data <- simDCM(range.mean.phi = c(0.8, 0.8), range.mean.gamma = c(0.3, 0.3),
  range.mean.p = c(0.6, 0.6)))

# Function arguments that lead to much structure (no zero arguments)
str(data <- simDCM(nspecies = 200, nsites = 267, nsurveys = 3, nyears = 25,
  mean.psi1 = 0.4, sig.lpsi1 = 3, mu.beta.lpsi1 = 1, sig.beta.lpsi1 = 3,
  range.mean.phi = c(0.5, 1), sig.lphi = 3, mu.beta.lphi = 1,
  sig.beta.lphi = 3, range.mean.gamma = c(0, 0.5),
  sig.lgamma = 3, mu.beta.lgamma = -1, sig.beta.lgamma = 3,
  range.mean.p = c(0.1, 0.9), sig.lp = 3, mu.beta.lp = 1, sig.beta.lp = 3,
  range.beta1.survey = c(-2, -0.5), range.beta2.survey = c(0, 2),
  trend.sd.site = c(0, 0), trend.sd.survey = c(0, 0), show.plot = TRUE))

# Not every occurring species will be detected
set.seed(1)
str(data <- simDCM(nspecies = 200, nsites = 20, nsurveys = 2, nyears = 10,
  mean.psi1 = 0.1, sig.lpsi1 = 5,
  range.mean.phi = c(0.3, 0.3), sig.lphi = 5,
  range.mean.gamma = c(0.1, 0.1), sig.lgamma = 5,
  range.mean.p = c(0.1, 0.1), sig.lp = 5) )

# Pull out data from species 5
ysp5 <- data$y[,,,5]

# Pull out data from year 1
yyr1 <- data$y[,,1,]
}
}
