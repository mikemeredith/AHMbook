\name{simDynoccSpatial}
\alias{simDynoccSpatial}
\encoding{UTF-8}
\title{
Simulate data under a dynamic occupancy model with spatial autocorrelation
}
\description{
Function to simulate detection/nondetection data under a general dynamic site-occupancy model with autocorrelation, including:

* Annual variation in the probabilities of patch persistence, colonization and detection is specified by the bounds of a uniform distribution.

* One covariate is allowed to affect a parameter: a site covariate for psi1, site-by-year covariates for phi and gamma, and an observational covariate for p. Covariates are generated internally from uniform(-2, 2) distributions.

* Additional heterogeneity among sites in persistence and colonization or both.

* Additional detection heterogeneity at the site or survey level, with the possibility of a temporal trend in this heterogeneity over the years. E.g., an annual trend in detection heterogeneity at the site or the survey level is specified by the first and second value, which correspond to the heterogeneity in the first and the last year. Hence, trend.sd.site = c(0, 1) will result in a linear trend in the magnitude of site heterogeneity in detection from 0 in the first year to 1 in the last year.

* A single, spatially structured covariate for habitat suitability may affect all parameters via coefficient beta.XAC; for a biologically reasonable way, choose coefficients with the same sign for all 4 (mediated by underlying density). That spatial covariate is simulated as a Gaussian random field with negative exponential correlation function with 'range parameter' theta.XAC.

* Autologistic effects (beta.Xautolog) in persistence and colonization probability can be chosen, which fits a logistic regression of these parameters on the proportion of occupied neighboring cells (in a queen's or 2nd order neighborhood) during the previous time step.

* Additional detection heterogeneity can be introduced at the site- or the individual survey level, with the possibility of a temporal trend in this heterogeneity. For instance, an annual trend in detection heterogeneity at the site or the survey level is specified by the value in the first and the last year. Hence, trend.sd.site = c(0, 1) will result in a linear trend in the magnitude of site heterogeneity in detection from 0 in the first year to 1 in the last year.

Previous versions used \pkg{RandomFields}, but that is not currently available on CRAN. If not available, \pkg{fields} will be used, but it cannot deal with large values of \code{side} and \code{theta.XAC}.

See also \code{\link{simDynoccSpatialData}}.
}
\usage{
simDynoccSpatial(side = 50, nyears = 10, nsurveys = 3,
      mean.psi1 = 0.4, beta.Xpsi1 = 0,
      range.phi = c(0.8, 0.8), beta.Xphi = 0,
      range.gamma = c(0.1, 0.1), beta.Xgamma = 0,
      range.p = c(0.4, 0.4), beta.Xp = 0,
      theta.XAC = 5000, beta.XAC = c(0, 0, 0, 0), beta.Xautolog = c(0, 0),
      trend.sd.site = c(0, 0), trend.sd.survey = c(0, 0),
      seed.XAC = NA, seed = NULL, show.plots = TRUE, ask.plot = TRUE, verbose = TRUE)
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
bounds of uniform distribution from which annual probability of colonization is randomly drawn.
}
  \item{beta.Xgamma}{
coefficient of environmental covariate in probability of colonization.
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
vector of coefficients of autologistic covariate in the following order: persistence (phi), colonization (gamma); the autocovariate is computed at every season as the proportion of occupied cells in a queen's neighborhood around each cell.
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
  \item{show.plots}{
if TRUE, summary plots are displayed.
}
  \item{ask.plot}{
If TRUE, pause between plots of results; set to FALSE if running simulations.
}
  \item{verbose}{
if TRUE, output will be written to the console.
}
}

\value{
A list with the values of the arguments input and the following additional elements:
  \item{M }{scalar, total number of pixels in the study area, side^2}
  \item{grid }{2-column matrix, x and y coordinates of each pixel}
  \item{amatrix }{M x M matrix, [i,j] = 1 if cells i and j are neighbors, 0 otherwise}
  \item{Xpsi1 }{side x side matrix, value of covariate affecting initial occupancy}
  \item{Xphi }{side x side x nyears array, value of covariate affecting persistence}
  \item{Xgamma }{side x side x nyears array, value of covariate affecting colonization}
 \item{Xp }{side x side x nsurveys x nyears array, value of covariate affecting detection}
  \item{XAC }{side x side matrix, the spatially correlated covariate}
 \item{Xauto }{side x side x nyears array, the autocovariate, the proportion of neighboring cells occupied}
 \item{Xautoobs }{side x side x nyears array, the observed autocovariate, the proportion of neighboring cells where the species was detected}
  \item{sd.site }{vector nyears, year-specific values of SD of Gaussian random site effects in p}
  \item{sd.survey }{vector nyears, year-specific values of SD of Gaussian random survey effects in p}
  \item{mean.phi }{vector nyears-1, mean persistence for each interval}
 \item{mean.gamma }{vector nyears-1, mean colonization for each interval}
 \item{mean.p }{vector nyears, mean detection probability for each year}
 \item{psi }{side x side x nyears array, probability of occupancy for each site and year}
  \item{mean.psi }{vector nyears, mean occupancy over cells for each year}
  \item{psi.app }{vector nyears, apparent occupancy, proportion of cells where species detected}
  \item{z }{side x side x years array, true occupancy state}
  \item{zobs }{side x side x years array, observed occupancy state}
  \item{nocc }{vector nyears, true number of occupied cells}
  \item{nocc.obs }{vector nyears, number of cells where the species was detected}
  \item{phi }{side x side x nyears-1 array, probability of persistence in each interval between years}
  \item{gamma }{side x side x nyears-1 array, probability of colonization in each interval between years}
  \item{p }{side x side x nsurveys x nyears array, probability of detection}
  \item{y }{side x side x nsurveys x nyears array, the observed detection history}
   \item{umf }{an unmarked data frame object with the simulated data}
}
\references{
Kéry, M. & Royle, J.A. (2021) \emph{Applied Hierarchical Modeling in Ecology} AHM2 - 9.

}
\author{
Marc Kéry & Andy Royle
}

\examples{
# Generate data with the default arguments and look at the structure
try(str(simDynoccSpatial()))  # Fails if RandomFields is not available

str(simDynoccSpatial(side=50, theta.XAC=1))

}
