\name{simOccSpatial}
\alias{simOccSpatial}
\encoding{UTF-8}
\title{
Simulate replicated detection/nondetection data under a spatial, static occupancy model
}
\description{
Simulates replicated detection/nondetection data under a spatial, static occupancy model for a semi-realistic landscape in a square of 50x50 km in the Bernese Oberland around Interlaken, Switzerland. Unit of the data simulation is a 1km2 quadrat, hence, there are 2500 units. It uses the data set \code{\link{BerneseOberland}}, which has covariates for elevation and forest cover.

For occupancy, the function allows you to specify a quadratic effect of elevation. Then, a Gaussian spatial random field (s) with negative exponential correlation function is simulated using \code{\link{simExpCorrRF}}; you can set the variance and the range (scale) parameter theta. Basically, the larger the value of theta.RF, the bigger are the 'islands' simulated in the random field. The abundance in each quadrat (i) is built up via the following linear predictor:

\code{psi[i] <- qlogis(beta0 + beta[1] * elev[i] + beta[2] * elev[i]^2 + s[i])} \cr
\code{z[i] ~ Bernoulli(psi[i])}

Replicated detection/nondetection data are simulated as usual under a Bernoulli observation model, and detection probability is allowed to vary by one site and one observational covariate: respectively quadrat forest cover (in the BerneseOberland data set), and wind-speed, which is invented data. Counts at each site (i) and for each occasion (j) are then produced according to the following model:

\code{p[i,j] <- plogis(alpha0 + alpha[1] * forest[i] + alpha[2] * wind[i,j])} \cr
\code{y[i,j] ~ Bernoulli(z[i] * p[i,j]) }

Finally, we assume that only a subset of the 2500 quadrats is surveyed. Hence, we allow you to choose the number of quadrats that are surveyed and these will then be randomly placed into the landscape. We then assume that counts will only be available for these surveyed quadrats, i.e., detection/non-detection data from all non-surveyed quadrats will be replaced with NA.
}
\usage{
simOccSpatial(nsurveys = 3, mean.psi = 0.6, beta = c(2, -2),
  mean.p = 0.4, alpha = c(-1, -1), sample.size = 500, variance.RF = 1, theta.RF = 10,
  seeds = c(10, 100), show.plots=TRUE)
}
\arguments{
  \item{nsurveys}{
number of surveys per quadrat.
}
  \item{mean.psi}{
probability of occupancy at a quadrat with mean values of all covariates; \code{beta0} is \code{qlogis(mean.psi)}.
}
  \item{beta}{
vector of length 2, the linear and quadratic coefficients for the effect of elevation on occupancy.
}
  \item{mean.p}{
Expected detection at the average value of all detection covariates (and ignoring all random effects): detection model
}
  \item{alpha}{
vector of length 2, the coefficients for the effect of forest and wind on detection.
}
  \item{sample.size}{
the number of quadrats surveyed.
}
  \item{variance.RF}{
variance of the random field.
}
  \item{theta.RF}{
parameter governing spatial correlation (=1/phi); the larger the value of theta.RF, the bigger are the 'islands' simulated in the random field.
}
  \item{seeds}{
vector of length 2; random seeds used for the random field and the selection of quadrats surveyed respectively.
}
  \item{show.plots}{
if TRUE, plots of the data will be displayed; set to FALSE if you are running many simulations.
}
}

\value{
A list with the arguments input and the following additional elements:
 \item{xcoord, ycoord }{The x and y coordinates from the \code{BerneseOberland} data set.}
 \item{elevation }{The elevation covariate from the \code{BerneseOberland} data set.}
 \item{forest }{The forest cover covariate from the \code{BerneseOberland} data set.}
 \item{elevationS }{The elevation covariate standardized to mean 0, SD 1.}
 \item{forestS }{The forest cover covariate standardized to mean 0, SD 1.}
 \item{wind }{The wind covariate, generated internally.}
 \item{field }{The spatially-autocorrelated covariate, generated internally.}
 \item{alpha0 }{Intercept for the detection model, equal to \code{qlogis(mean.p)}.}
 \item{beta0 }{Intercept for the abundance model, equal to \code{qlogis(mean.psi)}.}
 \item{z }{True occupancy in each quadrat.}
 \item{trueNocc }{True number of occupied sites.}
 \item{obsNocc }{Number of sites where the species was detected.}
 \item{true_psi_fs }{ }
 \item{obs_psi_fs }{ }
 \item{p }{probability of detection for each survey.}
 \item{y }{detected/non-detected data for each quadrat.}
 \item{surveyed.sites }{the IDs of the quadrats surveyed.}
 \item{yobs }{detected/non-detected data for each quadrat surveyed, NA for quadrats not surveyed.}
}
\references{
Kéry, M. & Royle, J.A. (2020) \emph{Applied Hierarchical Modeling in Ecology} AHM2 - 9.2.

}
\author{
Marc Kéry & Andy Royle
}
\examples{
# Generate data with the default arguments and look at the structure:
str(dat <- simOccSpatial())
str(dat <- simOccSpatial(show.plots=FALSE))
\donttest{
# More surveys
str(dat<- simOccSpatial(nsurveys = 10))

# Minimal number of surveys is 1
str(dat<- simOccSpatial(nsurveys = 1))

# A truly ubiquitous species
str(dat <- simOccSpatial(mean.psi = 1))

# Only negative linear effect of elevation
str(dat <- simOccSpatial(beta = c(2, 0)))

# No effect of elevation at all (see effects of spatial field now clearly)
str(dat <- simOccSpatial(beta = c(0, 0)))

# Perfect detection (p = 1)
str(dat <- simOccSpatial(mean.p = 1))

# No effect in detection of forest cover
str(dat <- simOccSpatial(alpha = c(0, -1)))

# No effect in detection of wind speed (see neatly forest effect now)
str(dat <- simOccSpatial(alpha = c(-1, 0)))

# Sample only 100 quadrats
str(dat <- simOccSpatial(sample.size = 100))

# Sample all 2500 quadrats
str(dat <- simOccSpatial(sample.size = 2500))

# Larger variance of the multivariate Gaussian Random variable in the random field
#  (this will increase the effect of the field on occupancy and detection)
str(dat <- simOccSpatial(variance.RF = 10))

# No spatial autocorrelation (Variant 1: set variance to 0)
str(dat <- simOccSpatial(variance.RF = 0))

# No spatial autocorrelation (Variant 2: set theta very close to 0,
#  but not quite 0, otherwise function breaks)
str(dat <- simOccSpatial(theta.RF = 0.0001))

# Larger value of theta.RF gives larger 'islands'
str(dat <- simOccSpatial(theta.RF = 100))

# Even larger value of theta.RF gives even larger 'islands'
str(dat <- simOccSpatial(theta.RF = 10000))
}
}
