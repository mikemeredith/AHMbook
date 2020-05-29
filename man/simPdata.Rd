\name{sim.pdata}
\alias{sim.pdata}
\encoding{UTF-8}
\title{
Simulate non-hierarchical point transect (= point count) data
}
\description{
Function simulates coordinates of individuals on a square with a count location at the center point.
}
\usage{
sim.pdata(N = 1000, sigma = 1, B = 3, keep.all = FALSE, show.plot = TRUE)
}
\arguments{
  \item{N}{
total population size in the square
}
  \item{sigma}{
scale of half-normal detection function
}
  \item{B}{
circle radius; the data are simulated over a square of side 2 * B, but individuals outside the circle of radius B are not detected.
}
  \item{keep.all}{
if TRUE, the data for all individuals are returned; if FALSE, only for individuals detected.
}
  \item{show.plot}{
choose whether to show plots or not. Set to FALSE when using function in simulations.
}
}

\value{
A list with the following elements:
 \item{N }{total population size in the square}
 \item{sigma }{scale of half-normal detection function}
 \item{B }{circle radius}
 \item{u1, u2 }{the x and y coordinates of each of the individuals}
 \item{d }{the distance of each individual from the center of the circle}
 \item{y }{a 0/1 indicator showing whether each individual is detected or not, a vector of length N}
 \item{N.real }{the realized number of individuals within the circle of radius B}

}
\references{
Kéry, M. & Royle, J.A. (2016) \emph{Applied Hierarchical Modeling in Ecology} AHM1 - 8.2.5.1.
}
\author{
Marc Kéry & Andy Royle
}

\examples{
# Simulate a data set with the default arguments and look at the structure of the output:
tmp <- sim.pdata()
str(tmp)
}
