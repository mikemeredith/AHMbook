\name{simHDSpoint}
\alias{simHDSpoint}
\encoding{UTF-8}

\title{
Simulate data under hierarchical distance sampling point transect protocol
}
\description{
The function simulates hierarchical distance sampling (HDS) data under a point transect protocol. At each site, it works with a circle of radius \code{B}.

The state process is simulated by first drawing a covariate value, "habitat", for each site from a Normal(0, 1) distribution. This is used in a log-linear regression with arguments \code{mean.density} and \code{beta.density} to calculate the expected density of animals per hectare. The expected number of animals in the circle is calculated from the area of the circle and the density, and numbers are drawn from a Poisson distribution for each site.

Animals are assumed to be distributed randomly over the circle, and distances from the point are generated.

A detection covariate, "wind", for each site is drawn from a Uniform(-2, 2) distribution. This is used in a log-linear regression with arguments \code{mean.sigma} and \code{beta.sigma} to calculate the scale parameter, sigma, of the half-normal detection function. Detections are simulated as Bernoulli trials with probability of success decreasing with distance from the point.

This is a simplified (and faster) version of the function simHDS with type="point", but works with densities and with measurements in hectares and meters.
}
\usage{
simHDSpoint(nsites = 1000, mean.density = 1, beta.density  = 1,
    mean.sigma = 20, beta.sigma = -5, B = 60, discard0=FALSE, show.plots=TRUE)
}
\arguments{
   \item{nsites}{
Number of sites (spatial replication)
}
  \item{mean.density}{
the expected value of density (animals per hectare) when the habitat covariate = 0; the intercept of the log-linear regression for density is log(mean.density).
}
  \item{beta.density}{
the slope of the log-linear regression for density on a habitat covariate.
}
  \item{mean.sigma}{
the expected value of the scale parameter of the half-normal detection function when the wind speed = 0; the intercept of the log-linear regression for sigma is log(mean.sigma).
}
  \item{beta.sigma}{
the slope of log-linear regression of scale parameter of the half-normal detection function on wind speed.
}
  \item{B}{
the circle radius in meters.
}
  \item{discard0}{
If TRUE, subset to sites at which individuals were captured. You may or may not want to do this depending on how the model is formulated so be careful.
}
  \item{show.plots}{
choose whether to show plots or not. Set to FALSE when using function in simulations.
}
}

\value{
A list with the values of the arguments entered and the following additional elements:

\item{data }{simulated distance sampling data: a matrix with a row for each individual detected and 2 columns: site ID and distance from the point in meters. If \code{discard0 = FALSE}, the default, sites with no detections will appear in the matrix with NAs in column 2.}
\item{counts }{the number of individuals detected at each site}
\item{habitat }{simulated habitat covariate}
\item{wind }{simulated detection covariate}
\item{N }{simulated number of individuals within the circle sampled}

}
\references{
Kéry, M. & Royle, J.A. (2016) \emph{Applied Hierarchical Modeling in Ecology} AHM1 - 8.5.1.
}
\author{
Marc Kéry, Andy Royle, Mike Meredith
}

\examples{
# Simulate a data set with the default arguments and look at the structure of the output
set.seed(123)
tmp <- simHDSpoint()
str(tmp)
head(tmp$data, 10)
head(tmp$counts, 10)

# Without 'wind', plots clearly show effect of distance from point
tmp <- simHDSpoint(beta.sigma=0)

tmp <- simHDSpoint(discard0=TRUE)
str(tmp)
head(tmp$data, 10) # some sites missing, no NAs.
head(tmp$counts)
}
