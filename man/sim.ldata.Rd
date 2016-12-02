\name{sim.ldata}
\alias{sim.ldata}
\encoding{UTF-8}

\title{
Simulation of distance sampling data.
}
\description{
Simulates non-hierarchical line transect data under CDS. It subjects N individuals to sampling, and then retains the value of distance from transect only for individuals that are captured.
}
\usage{
sim.ldata(N = 200, sigma = 30)
}
\arguments{
  \item{N}{
number of individuals along transect with distance u(-100, 100)
}
  \item{sigma}{
scale parameter of half-normal detection function
}
}

\value{
A listwith the following elements:
 \item{N }{the number of individuals along the transect}
 \item{sigma }{scale parameter of half-normal detection function}
 \item{xall }{distance from the transect line for all N individuals}
 \item{x }{absolute distance from the transect line for those individuals detected}
}
\references{
Kéry & Royle (2016) \emph{Applied Hierarachical Modeling in Ecology} Section 8.2.3
}
\author{
Kéry & Royle
}

\examples{
# Simulate a data set with the default arguments and look at the structure of the output:
tmp <- sim.ldata()
str(tmp)
}
