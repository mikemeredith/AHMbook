\name{simDSM}
\alias{simDSM}
\encoding{UTF-8}

\title{
Simulate line transect data for density surface modeling
}
\description{
The function generates a population represented as a binomial point pattern in a heterogeneous landscape with density a function of the covariate Habitat. Data for multiple line transect surveys using a wiggly transect are then simulated, and the pixel IDs for the activity centers of detected individuals returned.

To recreate the data sets used in the book with R 3.6.0 or later, include \code{sample.kind="Rounding"} in the call to \code{set.seed}. This should only be used for reproduction of old results.
}
\usage{
simDSM(X, Ntotal = 400, sigma = 0.65, beta1 = 1.0,
    nsurveys = 2, xlim = c(-0.5, 3.5), ylim = c(-0.5, 4.5), show.plots = TRUE)
}

\arguments{
  \item{X}{
a 2-column matrix with coordinates of regularly spaced points along the transect line; see Examples.
}
  \item{Ntotal}{
the true total number of individuals in the study area.
}
  \item{sigma}{
scale parameter for the half-normal detection function.
}
  \item{beta1}{
coefficient for the relationship between the Habitat covariate and population density.
}
  \item{nsurveys}{
the number of replicate surveys along the transect.
}
  \item{xlim, ylim}{
the extent of the (rectangular) study area
}
  \item{show.plots}{
if TRUE, summary plots are displayed.
}
}

\value{
A list with the values of the input arguments and the following additional elements:
 \item{Habitat }{a vector for the habitat covariate for each pixel}
 \item{Habgrid }{a 2-column matrix with the coordinates of center of each pixel}
 \item{nPix }{the number of pixels in the study area}
 \item{N }{true number of activity centers in each pixel}
 \item{U }{a 2-column matrix with the locations of ACs for all individuals in the population}
 \item{Ucap }{a 2-column matrix with the locations of ACs for individuals detected at least once}
 \item{nind }{the number of individuals detected at least once}
 \item{pixel }{a nind x nsurvey matrix with the pixel ID for the activity center or NA if the individual was not detected on the survey}
}
\references{
Kéry, M. & Royle, J.A. (2021) \emph{Applied Hierarchical Modeling in Ecology} AHM2 - 11.10.
}
\author{
Marc Kéry, Andy Royle & Mike Meredith
}

\examples{
# Run the function with default values and look at the output
library(AHMbook)
data(wigglyLine)
points <- sp::SpatialPoints( wigglyLine )
sLine <- sp::Line(points)
regpoints <- sp::spsample(sLine, 100, type = "regular")
str(simDSM(X = regpoints@coords))

# Generate the data set used in AHM2 11.10
RNGversion("3.5.3")
set.seed(2027, kind = "Mersenne-Twister")
tmp <- simDSM(X = regpoints@coords) # Produces Fig 11.15 in the book
}
