\name{sim.fn}
\alias{sim.fn}
\encoding{UTF-8}
\title{
Simulate a Poisson point process
}
\description{
Simulates animal or plant locations in space according to a homogenous Poisson process. This process is characterized by the intensity, which is the average number of points per (very small)unit area. The resulting point pattern is then discretized to obtain abundance data and  presence/absence (or occurrence) data. The discretization of space is achieved by choosing the cell size. It is used in AHM1 Section 1.1 to help to understand the relationship between point patterns, abundance data and occurrence data (also called presence/absence or distribution data). For a similar and somewhat more sexy version of this function, see \code{\link{simPPe}}.
}
\usage{
sim.fn(quad.size = 10, cell.size = 1, intensity = 1, show.plot = TRUE)
}
\arguments{
  \item{quad.size}{
The length of each side of the quadrat (in arbitrary units)
}
  \item{cell.size}{
The length of each side of the cells into which the quadrat is divided. The ratio of quad.size to cell.size must be an integer.
}
  \item{intensity}{
The average number of points (animals or plants) per unit area.
}
  \item{show.plot}{
If TRUE, the results are plotted. Set to FALSE when running simulations.
}
}
\value{
A list with the values of the arguments and the following additional elements:
\item{exp.N }{Expected population size in quadrat}
\item{breaks }{boundaries of grid cells}
\item{n.cell }{Number of cells in the quadrat}
\item{mid.pt }{Cell mid-points}
\item{M }{ Realized population size in quadrat}
\item{u1 }{ x coordinate of each individual}
\item{u2 }{ y coordinate of each individual}
\item{N }{The number of individuals in each cell (a vector of length n.cell)}
\item{z }{Presence/absence (1/0) in each cell (a vector of length n.cell)}
\item{psi }{Proportion of cells occupied, ie, the species is present.}

}
\references{
Kéry, M. & Royle, J.A. (2016) \emph{Applied Hierarchical Modeling in Ecology} AHM1 - 1.1.
}
\author{
Marc Kéry and Andy Royle
}

\examples{
# Generate a simulated data set with default arguments and look at the structure:
tmp <- sim.fn()
str(tmp)

# Effect of grain size of study on abundance and occupancy (intensity constant)
tmp <- sim.fn(quad.size = 10, cell.size = 1, intensity = 0.5)
tmp <- sim.fn(quad.size = 10, cell.size = 2, intensity = 0.5)
tmp <- sim.fn(quad.size = 10, cell.size = 5, intensity = 0.5)
tmp <- sim.fn(quad.size = 10, cell.size = 10, intensity = 0.5)

# Effect of intensity of point pattern (intensity) on abundance and occupancy
tmp <- sim.fn(intensity = 0.1) # choose default quad.size = 10, cell.size = 1
tmp <- sim.fn(intensity = 1)
tmp <- sim.fn(intensity = 5)
tmp <- sim.fn(intensity = 10)
}
