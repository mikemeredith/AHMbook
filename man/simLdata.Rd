\name{sim.ldata}
\alias{sim.ldata}
\encoding{UTF-8}

\title{
Simulation of distance sampling data.
}
\description{
Simulates non-hierarchical line transect data under conventional distance sampling (CDS). It subjects N individuals to sampling, and then retains the value of distance from transect only for individuals that are captured.
}
\usage{
sim.ldata(N = 200, sigma = 30, show.plot = TRUE)
}
\arguments{
  \item{N}{
number of individuals along transect with distance Uniform(-100, 100)
}
  \item{sigma}{
scale parameter of half-normal detection function
}
  \item{show.plot}{
choose whether to show plots or not. Set to FALSE when using function in simulations.
}
}

\value{
A list with the following elements:
 \item{N }{the number of individuals along the transect}
 \item{sigma }{scale parameter of half-normal detection function}
 \item{xall }{distance from the transect line for all N individuals}
 \item{x }{absolute distance from the transect line for those individuals detected}
}
\references{
Kéry, M. & Royle, J.A. (2016) \emph{Applied Hierarchical Modeling in Ecology} AHM1 - 8.2.3.
}
\author{
Marc Kéry & Andy Royle
}

\examples{
# Simulate a data set with the default arguments and look at the structure of the output:
tmp <- sim.ldata()
str(tmp)
}
