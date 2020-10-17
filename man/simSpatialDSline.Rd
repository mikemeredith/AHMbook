\name{simSpatialDSline}
\alias{simSpatialDSline}
\encoding{UTF-8}

\title{
Simulate data for a standardized line transect
}
\description{
This simulates line transect distance sampling data with a spatial distribution of objects in a heterogeneous landscape where the parameter beta controls the effect of habitat. Habitat is simulated according to a Gaussian random field model defined within the function. Uses a half normal detection model (if perp = TRUE) or a Gaussian hazard model (perp = FALSE).

To recreate the data sets used in the book with R 3.6.0 or later, include \code{sample.kind="Rounding"} in the call to \code{set.seed}. This should only be used for reproduction of old results.
}
\usage{
simSpatialDSline(N=1000, beta = 1, sigma=0.25, alpha0 = -2, W=1/2, L = 4,
    perp=FALSE, show.plots=TRUE)
}

\arguments{
  \item{N}{
total population size in the rectangle
}
  \item{beta}{
coefficient of spatial covariate x for the density model.
}
  \item{sigma}{
scale of half-normal detection function
}
  \item{alpha0}{
intercept of the hazard function.
}
\item{W}{
half-width of the rectangle, which extends W each side of the transect line.
}
  \item{L}{
length of the transect.
}
  \item{perp}{
if TRUE, data are simulated for a traditional distance sampling model with perpendicular distances; if FALSE (the default) a model with 'forward distance' data, ie, the distance from the observer to the animal on first detection.
}
  \item{show.plots}{
if TRUE, summary plots are displayed.
}
}

\value{
A list with the values of the input arguments and the following additional elements:
 \item{delta }{the distance between pixel centers (spatial resolution of the raster}
 \item{grid }{2-column matrix with x/y coordinates of all pixels}
 \item{Habitat }{value of habitat covariate for each pixel}
 \item{Habraster }{a Raster object with the habitat covariate}
 \item{u1, u2 }{x and y coordinates for all the animals in the population}
 \item{traps }{2-column matrix of trap locations}
If \code{perp = TRUE} we have
 \item{data }{a 2-column matrix with x and y coordinates of each animal captured.}
 \item{pixel }{pixel ID for each animal captured.}
and if \code{perp = FALSE} we have
 \item{data }{a matrix with rows for each animal captured and columns for trap of first capture, distance from trap to animal, and x and y coordinates of the animal.}
 \item{pbar }{probability that each animal is the population is captured at least once}
 \item{pixel }{pixel ID for each animal captured.}
}
\references{

Kéry, M. & Royle, J.A. (2021) \emph{Applied Hierarchical Modeling in Ecology} AHM2 - 11.
}
\author{
Marc Kéry & Andy Royle
}

\examples{
# Run the function with default values and look at the output
str(tmp <- simSpatialDSline(), 1)  # use str(., max.level=1) to limit the amount of output.
}
