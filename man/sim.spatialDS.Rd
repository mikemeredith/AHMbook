\name{sim.spatialDS}
\alias{sim.spatialDS}
\encoding{UTF-8}

\title{
Simulates data under spatial hierarchical distance sampling model
}
\description{
Function simulates coordinates of individuals on a square
}
\usage{
sim.spatialDS(N = 1000, beta = 1, sigma = 1, keep.all = FALSE, B = 3,
  model=c("halfnorm", "hazard"), show.plot=TRUE)
}
\arguments{
  \item{N}{
total population size in the square
}
  \item{beta}{
coefficient of SOMETHING on spatial covariate x
}
  \item{sigma}{
scale parameter of detection function
}
  \item{keep.all}{
if TRUE, the data for all individuals are returned; if FALSE, only for individuals detected.
}
  \item{B}{
circle radius; the data are simulated over a square of side 2 * B %, but individuals outside the circle of radius B are not detected. % In this case they ARE detected.
}
  \item{model}{
The detection function used, either "halfnorm" or "hazard".
}
  \item{show.plot}{
choose whether to show plots or not. Set to FALSE when using function in simulations.
}
}
\value{
A list with the values of the input arguments and the following additional elements:
  \item{u1 }{x coordinate of each animal}
  \item{u2 }{y coordinate of each animal}
  \item{d }{distance of each animal from the centre of the circle}
  \item{y }{indicator of detection of each animal, a vector of length N}
  \item{N.real }{the number of animals inside the circle of radius B}
  \item{Habitat }{Value of the spatially  correlated habitat covariate, a 900 x 1 matrix}
  \item{grid }{Coordinates of the centres of each pixel, a 900 x 2 matrix}
}
\references{
Kéry & Royle (2016) \emph{Applied Hierarachical Modeling in Ecology} Section 9.8.3

}
\author{
Kéry & Royle
}

\examples{
# Generate data with the default arguments and look at the structure:
tmp <- sim.spatialDS()
str(tmp)
}
