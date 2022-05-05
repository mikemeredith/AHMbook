\name{simNmixSpatial}
\alias{simNmixSpatial}
\encoding{UTF-8}
\title{
Simulate replicated counts under a spatial, static binomial N-mixture model
}
\description{
Simulates replicated counts under a spatial, static binomial N-mixture model for a semi-realistic landscape in a square of 50x50 km in the Bernese Oberland around Interlaken, Switzerland. Unit of the data simulation is a 1km2 quadrat, hence, there are 2500 units. It uses the data set \code{\link{BerneseOberland}}, which has covariates for elevation and forest cover.

For abundance, the function allows you to specify a quadratic effect of elevation. Then, a Gaussian spatial random field (s) with negative exponential correlation function is simulated using \code{\link{simExpCorrRF}}; you can set the variance and the range (scale) parameter theta. Basically, the larger the value of theta.RF, the bigger are the 'islands' simulated in the random field. The abundance in each quadrat (i) is built up via the following linear predictor:

\code{lambda[i] <- exp(beta0 + beta[1] * elev[i] + beta[2] * elev[i]^2 + s[i])} \cr
\code{N[i] ~ Poisson(lambda[i])}

Replicated counts are simulated as usual under a binomial observation model, and detection probability is allowed to vary by one site and one observational covariate: respectively quadrat forest cover (in the BerneseOberland data set), and wind-speed, which is invented data. Counts at each site (i) and for each occasion (j) are then produced according to the following model:

\code{  p[i,j] <- plogis(alpha0 + alpha[1] * forest[i] + alpha[2] * wind[i,j])} \cr
\code{C[i,j] ~ Binomial(N[i], p[i,j]) }

Finally, we assume that only a subset of the 2500 quadrats is surveyed. Hence, we allow you to choose the number of quadrats that are surveyed and these will then be randomly placed into the landscape. We then assume that counts will only be available for these surveyed quadrats, i.e., counts from all non-surveyed quadrats will be NA.

To recreate the data sets used in the book with R 3.6.0 or later, call \code{RNGversion="3.5.3"} before the call to \code{simNmixSpatial}. This should only be used for reproduction of old results.
}
\usage{
simNmixSpatial(nsurveys = 3, mean.lambda = exp(2), beta = c(2, -2),
  mean.p = 0.5, alpha = c(-1, -1), sample.size = 500, variance.RF = 1, theta.RF = 10,
  seeds = c(10, 100), truncN = 6, show.plots=TRUE, verbose = TRUE)
}
\arguments{
  \item{nsurveys}{
number of surveys per quadrat.
}
  \item{mean.lambda}{
expected number of individuals at a quadrat with mean values of all covariates; \code{beta0} is \code{log(mean.lambda)}.
}
  \item{beta}{
vector of length 2, the linear and quadratic coefficients for the effect of elevation on abundance.
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
  \item{truncN}{
used to limit the range of values of N displayed in some plots.
}
  \item{show.plots}{
 if TRUE, plots of the data will be displayed; set to FALSE if you are running simulations.
}
  \item{verbose}{
if TRUE, output will be written to the console.
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
 \item{beta0 }{Intercept for the abundance model, equal to \code{log(mean.lambda)}.}
 \item{lam }{Expected abundance in each quadrat.}
 \item{N }{Number of individuals at each site.}
 \item{Ntotal }{Total number of individuals in the study area.}
 \item{p }{probability of detection for each survey.}
 \item{y }{number of individuals detected in each quadrat.}
 \item{surveyed.sites }{the IDs of the quadrats surveyed.}
 \item{yobs }{number of individuals detected in each quadrat surveyed, NA for quadrats not surveyed.}
}
\references{
Kéry, M. & Royle, J.A. (2021) \emph{Applied Hierarchical Modeling in Ecology} AHM2 - 11.

}
\author{
Marc Kéry & Andy Royle
}
\examples{
# Generate data with the default arguments and look at the structure:
str(dat <- simNmixSpatial())
str(dat <- simNmixSpatial(show.plots=FALSE))
\donttest{
# More surveys
str(dat<- simNmixSpatial(nsurveys = 10))

# Minimal number of surveys is 1
str(dat<- simNmixSpatial(nsurveys = 1))

# Much more common species
str(dat<- simNmixSpatial(mean.lambda = exp(4)))

# Only negative linear effect of elevation
str(dat<- simNmixSpatial(beta = c(-2, 0)))

# No effect of elevation at all
str(dat<- simNmixSpatial(beta = c(0, 0)))

# Perfect detection (p = 1)
str(dat<- simNmixSpatial(mean.p = 1))

# No effect of forest cover on detection
str(dat<- simNmixSpatial(alpha = c(0, -1)))

# No effect of wind speed on detection
str(dat<- simNmixSpatial(alpha = c(-1, 0)))

# Sample only 100 quadrats
str(dat<- simNmixSpatial(sample.size = 100))

# Sample all 2500 quadrats
str(dat<- simNmixSpatial(sample.size = 2500))

# Larger variance of the multivariate Gaussian Random variable in the random field
#  (this will increase the effect of the field on abundance and counts)
str(dat<- simNmixSpatial(variance.RF = 10))

# No spatial autocorrelation (Variant 1: set variance to 0)
str(dat<- simNmixSpatial(variance.RF = 0))

# No spatial autocorrelation (Variant 2: set theta very close to 0,
#  but not quite 0, otherwise function breaks)
str(dat<- simNmixSpatial(theta.RF = 0.0001))

# Larger value of theta.RF gives larger 'islands'
try(str(dat<- simNmixSpatial(theta.RF = 100))) # works with RandomFields

# Truncate abundance in final plots to presence/absence
str(dat<- simNmixSpatial(truncN = 0.5))

# Essentially do not truncate abundance in final plots
str(dat<- simNmixSpatial(truncN = 70))
}
}
