\name{sim.spatialHDS}
\alias{sim.spatialHDS}
\encoding{UTF-8}
\title{
Simulates data for a hierarchical spatial distance sampling model
}
\description{
Generates data for distance sampling from spatially-replicated point transects, with density dependent on a spatially correlated habitat covariate. For each point count, the procedure is:

1. Simulate the habitat covariate over a grid of pixels covering a square.

2. Distribute the population of individuals over the square with probability of location in a pixel related to the covariate.

3. Decide which individuals are detected using a distance sampling model with an observer at the center of the square, with a half normal detection function. (Note that individuals outside the circle of radius \code{B} can be detected.)

The locations and detection status of individuals at all sites are collated and returned, except for individuals at sites when none are detected.

To recreate the data sets used in the book with R 3.6.0 or later, include \code{sample.kind="Rounding"} in the call to \code{set.seed}. This should only be used for reproduction of old results.
}
\usage{
sim.spatialHDS(lam0 = 4, sigma = 1.5, B = 3, nsites = 100,
  beta1 = 1, npix = 20, show.plots = 3)
}
\arguments{
  \item{lam0}{
the expected number of individuals in the square of side = 2*B if the habitat covariate had no effect (ie, \code{beta1 = 0}).
}
  \item{sigma}{
scale parameter of the half-normal detection function
}
  \item{B}{
distance from the observer to the side of the square. This is usually set so that the probability of detection of individuals outside the square is negligable, eg, \code{B = 3*sigma}.
}
  \item{nsites}{
number of sites
}
  \item{beta1}{
the size of the effect of the habitat covariate on the number of individuals in a pixel.
}
  \item{npix}{
the number of pixels along each dimension of the square: the entire grid has \code{npix x npix} pixels.
}
  \item{show.plots}{
the number of sites for which plots should be displayed.
}
}
\value{
A list with the following components:
 \item{data }{a matrix with columns for siteID, the coordinates of each individual (\code{u1, u2}), distance from the center of the square (\code{d}) and detection status, 0/1 (\code{y}), and rows for each individual, except that individuals at sites where no individuals were detected are NOT included, and those sites are represented by NAs}
 \item{B }{the radius of the circle; the data are simulated over a square of side 2*B}
 \item{Habitat }{a matrix with the values of the Habitat covariate for each pixel at each site}
 \item{grid }{a data frame with the coordinates of each pixel (same for all sites)}
 \item{N }{the realized number of individuals in the square at each site}
 \item{nsites }{the number sites}
}
\references{
Kéry, M. & Royle, J.A. (2016) \emph{Applied Hierarchical Modeling in Ecology} AHM1 - 9.8.5.

}
\author{
Marc Kéry & Andy Royle
}

\examples{
# Generate data with the default arguments and look at the structure:
tmp <- sim.spatialHDS()
str(tmp)
}
