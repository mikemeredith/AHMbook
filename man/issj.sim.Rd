\name{issj.sim}
\alias{issj.sim}
\encoding{UTF-8}
\title{
Simulate open distance sampling data for the Island Scrub Jays
}
\description{
Function to simulate open distance sampling data for the Island Scrub Jays, based on Sollmann et al (2015).
}
\usage{
issj.sim(B, db, lam, sigma, phi, gamma, npoints, nyrs, nbsize = -1.02)
}
\arguments{
  \item{B}{
Radius of the circle sampled.
}
  \item{db}{
Distance categories; a vector of cut points from 0 to B inclusive.
}
  \item{lam}{
Expected abundance per point count, a npoints x 1 matrix.
}
  \item{sigma}{
Scale parameter of the half-normal detection function at each point, a npoints x 1 matrix.
}
  \item{phi}{
Survival probability
}
  \item{gamma}{
Recruitment rate
}
  \item{npoints}{
Number of points
}
  \item{nyrs}{
Number of years
}
  \item{nbsize}{
Size parameter for the negative binomial distribution used to generate individual counts per point count circle for year 1.
}
}
\value{
A list with the following elements:
 \item{NcList }{A list with one element per year, with distances of all animals from the point.}
 \item{detList }{A list with one element per year, a npoints x nbands matrix with the number of animals detected in each distance category for each point.}
 \item{N }{The number of animals at each point for each year, a npoints x nyears matrix.}
 \item{cell }{???.}
 \item{parms }{A vector with the input parameters, \code{c(lam,phi, gamma, sigma)}}
 \item{y }{nsites x nyears matrix of total detections}
 \item{dclass }{a vector with the distance class for each animal detected}
 \item{site }{a corresponding vector with the site for each animal detected}
}
\references{
Kéry & Royle (2016) \emph{Applied Hierarachical Modeling in Ecology} Section 9.7.1

Sollmann, R., Gardner, B., Chandler, R.B., Royle, J.A., Sillett, T.S. (2015) An open population hierarchical distance sampling model. \emph{Ecology} 96, 325-331.
}
\author{
Based on Sollmann et al (2015)
}

\examples{
## TO DO. Need values for lam and sigma to be able to run it.
}
