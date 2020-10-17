\name{simCJS}
\alias{simCJS}
\encoding{UTF-8}

\title{
Simulate individual capture-histories under a Cormack-Jolly-Seber (CJS) survival model
}
\description{
Function generates individual capture-histories under a CJS model with possibly time-dependent parameters. The number of values for interval-specific survival (phi) and time-specific detection (p) must be ensured to be equal to the number of occasions (n.occ) minus 1.
}
\usage{
simCJS(n.occ = 6, n.marked = 20, phi = 0.7, p = 0.4, show.plot = TRUE)
}
\arguments{
  \item{n.occ}{
The number of occasions (eg, years).
}
  \item{n.marked}{
The number of individuals marked per year except for the last year; scalar, or a vector of length (n.occ-1).
}
  \item{phi}{
The apparent survival probability between years; scalar, or a vector of length (n.occ-1).
}
  \item{p}{
The recapture probability for all except the first year; scalar, or a vector of length (n.occ-1).
}
  \item{show.plot}{
Choose whether to show plots or not. Set to FALSE when using function in simulations.
}
}
\value{
A list with the values of the arguments entered and the following additional elements:

\item{z }{n.ind x n.occ matrix, 1 if the individual is alive and in the study area, 0 otherwise.}
\item{ch }{n.ind x n.occ matrix, 1 if the individual is captured, 0 otherwise.}
\item{f }{n.ind vector, occasion the individual was first captured, and marked.}
\item{n.ind }{scalar, total number of individuals marked.}
\item{n.alive }{n.occ vector, number of animals alive and in the study area each year.}
}
\references{
Kéry, M. & Schaub, M. (2012) \emph{Bayesian population analysis using WinBUGS - a hierarchical perspective}, Academic Press.

Kéry, M. & Royle, J.A. (2021) \emph{Applied Hierarchical Modeling in Ecology} AHM2 - 3.2.2.
}
\author{
Marc Kéry & Andy Royle, based on code written by Michael Schaub for Chapter 7 of Kéry & Schaub (2012).
}

\examples{
# Run with the default arguments and look at the structure of the output:
set.seed(123)
tmp <- simCJS()
str(tmp)
colSums(tmp$ch)
}
