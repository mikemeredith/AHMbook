\name{simDynoccSpatial}
\alias{simDynoccSpatial}
\encoding{UTF-8}
\title{
Simulate data under a dynamic occupancy model with spatial components
}
\description{
Function to simulate detection/nondetection data under a general dynamic site-occupancy model, including:

* Annual variation in the probabilities of patch persistence, colonization and detection is specified by the bounds of a uniform distribution.

* One covariate is allowed to affect a parameter: a site covariate for psi1, site-by-year covariates for phi and gamma, and an observational covariate for p. Covariates are generated internally from uniform(-2, 2) distributions.

* Additional heterogeneity among sites in persistence and colonization or both.

* Additional detection heterogeneity at the site or survey level, with the possibility of a temporal trend in this heterogeneity over the years. E.g., an annual trend in detection heterogeneity at the site or the survey level is specified by the first and second value, which correspond to the heterogeneity in the first and the last year. Hence, trend.sd.site = c(0, 1) will result in a linear trend in the magnitude of site heterogeneity in detection from 0 in the first year to 1 in the last year.

* A single, spatially structured covariate for habitat suitability may affect all parameters via coefficient beta.XAC (for a biologically reasonable way, choose coefficients with the same sign for all 4 (mediated by underlying density). That spatial covariate is simulated as a Gaussian random field with negative exponential correlation function with 'range parameter' theta.XAC.

* Autologistic effects (beta.Xautolog) in persistence and colonization probability can be chosen, which fits a logistic regression of these parameters on the proportion of occupied neighbouring cells (in a queen's or 2nd order neighbourhood) during the previous time step.

* Additional detection heterogeneity can be introduced at the site- or the individual survey level, with the possibility of a temporal trend in this heterogeneity. For instance, an annual trend in detection heterogeneity at the site or the survey level is specified by the value in the first and the last year. Hence, trend.sd.site = c(0, 1) will result in a linear trend in the magnitude of site heterogeneity in detection from 0 in the first year to 1 in the last year.
}
\usage{
simDynoccSpatial(side = 50, nyears = 10, nsurveys = 3,
      mean.psi1 = 0.4, beta.Xpsi1 = 0,
      range.phi = c(0.8, 0.8), beta.Xphi = 0,
      range.gamma = c(0.1, 0.1), beta.Xgamma = 0,
      range.p = c(0.4, 0.4), beta.Xp = 0,
      theta.XAC = 5000, beta.XAC = c(0, 0, 0, 0), beta.Xautolog = c(0, 0),
      trend.sd.site = c(0, 0), trend.sd.survey = c(0, 0),
      seed.XAC = NA, seed = NULL, ask.plot = TRUE)
}
\arguments{
  \item{side}{
side length of square simulation area; the number of sites, or cells, M = side^2.
}
  \item{nyears}{
Number of years (or 'seasons').
}
  \item{nsurveys}{
Number of replicate surveys (= occasions) within a year.
}
  \item{mean.psi1}{
average occupancy probability in first year.
}
  \item{beta.Xpsi1}{
coefficient of environmental covariate in probability of initial occupancy.
}
  \item{range.phi}{
bounds of uniform distribution from which annual probability of persistence is randomly drawn.
}
  \item{beta.Xphi}{
coefficient of environmental covariate in probability of persistence.
}
  \item{range.gamma}{
bounds of uniform distribution from which annual probability of colonisation is randomly drawn.
}
  \item{beta.Xgamma}{
coefficient of environmental covariate in probability of colonisation.
}
  \item{range.p}{
bounds of uniform distribution from which probability of detection is randomly drawn.
}
  \item{beta.Xp}{
coefficient of environmental covariate in probability of detection.
}
  \item{theta.XAC}{
'range parameter' of a covariate with exponential spatial correlation (i.e., a Gaussian random field is used as an environmental covariate); must be > 0.
}
 \item{beta.XAC}{
vector of coefficients of that field for the 4 model parameters: psi1, phi, gamma, and p (in that order).
}
 \item{beta.Xautolog}{
vector of coefficients of autologistic covariate in the following order: persistence (phi), colonization (gamma); the autocovariate is computed at every season as the proportion of occupied cells in a queen's neighbourhood around each cell.
}
  \item{trend.sd.site}{
initial and final values of sd of normal distribution to model logit-normal noise in p at the site level; a linear trend is assumed over time; if the two values are the same, a constant value is assumed.
}
  \item{trend.sd.survey}{
initial and final values of sd of normal distribution to model logit-normal noise in p at the 'survey' level; if they are different, a linear trend is assumed over time.
}
  \item{seed.XAC}{
the seed to be used for simulation of the spatially structured covariate for habitat suitability.
}
  \item{seed}{
the seed to be used for simulation of values apart from the spatially structured covariate for habitat suitability; using the same value for seed.XAC with different values for seed allows generation of different data sets with the same habitat suitability covariate.
}
  \item{ask.plot}{
If TRUE, pause between plots of results; set to FALSE if running many simulations.
}
}

\value{
A list with the values of the arguments input and the following additional elements:
  \item{M }{scalar, total number of pixels in the study area, side^2}
  \item{grid }{2-column matrix, x and y coordinates of each pixel}
  \item{amatrix }{M x M matrix, [i,j] = 1 if cells i and j are neighbours, 0 otherwise}
  \item{Xpsi1 }{side x side matrix, value of covariate affecting initial occupancy}
  \item{Xphi }{side x side x nyears array, value of covariate affecting persistence}
  \item{Xgamma }{side x side x nyears array, value of covariate affecting colonisation}
 \item{Xp }{side x side x nsurveys x nyears array, value of covariate affecting detection}
  \item{XAC }{side x side matrix, the spatially correlated covariate}
 \item{Xauto }{side x side x nyears array, the autocovariate, the proportion of neighbouring cells occupied}
 \item{Xautoobs }{side x side x nyears array, the observed autocovariate, the proportion of neighbouring cells where the species was detected}
  \item{sd.site }{vector nyears, year-specific values of SD of Gaussian random site effects in p}
  \item{sd.survey }{vector nyears, year-specific values of SD of Gaussian random survey effects in p}
  \item{mean.phi }{vector nyears-1, mean persistence for each interval}
 \item{mean.gamma }{vector nyears-1, mean colonisation for each interval}
 \item{mean.p }{vector nyears, mean detection probability for each year}
 \item{psi }{side x side x nyears array, probability of occupancy for each site and year}
  \item{mean.psi }{vector nyears, mean occupancy over cells for each year}
  \item{psi.app }{vector nyears, apparent occupancy, proportion of cells where species detected}
  \item{z }{side x side x years array, true occupancy state}
  \item{zobs }{side x side x years array, observed occupancy state}
  \item{nocc }{vector nyears, true number of occupied cells}
  \item{nocc.obs }{vector nyears, number of cells where the species was detected}
  \item{phi }{side x side x nyears-1 array, probability of persistence in each interval between years}
  \item{gamma }{side x side x nyears-1 array, probability of colonisation in each interval between years}
  \item{p }{side x side x nsurveys x nyears array, probability of detection}
  \item{y }{side x side x nsurveys x nyears array, the observed detection history}
   \item{umf }{an unmarked data frame object with the simulated data}
}
\references{
Kéry & Royle (2019) \emph{Applied Hierarchical Modeling in Ecology} Chapter 9

}
\author{
Marc Kéry, 2014-2018
}

\examples{
# Generate data with the default arguments and look at the structure:
tmp <- simDynoccSpatial()
str(tmp)

\donttest{  # rather slow
# Large grid, seed.XAC = 88 gives a cool pattern
str(tmp <- simDynoccSpatial(side = 100, seed.XAC = 88))

# Longer series and no browsing through plots
str(tmp <- simDynoccSpatial(nyears = 100, ask.plot = FALSE))

# Dynamic system (high variance between years)
str(tmp <- simDynoccSpatial(range.phi = c(0, 0.8), range.gamma = c(0, 0.6)))

# Time-variation in all potentially time-varying parameters
str(tmp <- simDynoccSpatial(range.phi = c(0.5, 1), range.gamma = c(0, 0.5),
  range.p = c(0.1, 0.9)))

# Strong spatial autocorrelation via a spatially structured covariate
str(tmp <- simDynoccSpatial(beta.XAC = c(1, 1, 1, 1), beta.Xautolog = c(0, 0)))

# Strong spatial autocorrelation via a spatially structured covariate
str(tmp <- simDynoccSpatial(beta.XAC = c(1, 1, 1, 1), beta.Xautolog = c(2, 2)))

# Strong spatial autocorrelation via a spatially structured covariate
# but now effect of that spatial field only on colonization and persistence
str(tmp <- simDynoccSpatial(beta.XAC = c(0, 1, 1, 0), beta.Xautolog = c(0, 0)))

# Can it ever go extinct ???
str(tmp <- simDynoccSpatial(side = 20, nyears = 20, mean.psi1 = 0.01,
  range.phi = c(0.1, 0.1), range.gamma = c(0.1, 0.1), beta.Xautolog = c(0, 0)))

# Can it ever go extinct ???   -------> Yes
str(tmp <- simDynoccSpatial(side = 20, nyears = 10, mean.psi1 = 0.01,
  range.phi = c(0.3, 0.3), range.gamma = c(0, 0), beta.Xautolog = c(0, 0)), seed = 24)

# Evolution of number of occupied quads
apply(tmp$z, 3, sum)
# [1] 2 0 0 0 0 0 0 0 0 0
}
}
