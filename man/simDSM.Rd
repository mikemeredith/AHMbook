\name{simDSM}
\alias{simDSM}
\encoding{UTF-8}

\title{
Simulate line transect data for density surface modeling
}
\description{
The function generates a population represented as a point pattern in a heterogeneous landscape,
ie, a Poisson point pattern (PPP) with intensity a function of the covariate Habitat. Data for multiple line transect surveys using a wiggly transect are then simulated, and the pixel IDs for the activity centres of detected individuals returned.
}
\usage{
simDSM(X, N = 400, sigma = 0.65, beta1 = 1.0,
    nsurveys = 2, xlim = c(-0.5, 3.5), ylim = c(-0.5, 4.5), show.plot = TRUE)
}

\arguments{
  \item{X}{
a 2-column matrix with coordinates of _regularly_spaced_ points along the transect line; see Examples.
}
  \item{N}{
the true number of individuals in the study area.
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
  \item{show.plot}{
if TRUE, summary plots are displayed.
}
}

\value{
A list with the values of the input arguments and the following additional elements:
 \item{Habitat }{a vector for the habitat covariate for each pixel}
 \item{Habgrid }{a 2-column matrix with the coordinates of each pixel}
 \item{nPix }{the number of pixels in the study area}
 \item{nind }{the number of individuals detected at least once}
 \item{pixel }{a nind x nsurvey matrix with the pixel ID for the activity centre or NA if the individual was not detected on the survey}
}
\references{
Kéry, M. & Royle, J.A. (2020) \emph{Applied Hierarchical Modeling in Ecology} AHM2 - 11.10.
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
tmp <- simDSM(X = regpoints@coords)
str(tmp)
}
