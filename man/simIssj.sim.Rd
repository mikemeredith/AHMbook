\name{issj.sim}
\alias{issj.sim}
\encoding{UTF-8}
\title{
Simulate open distance sampling data for the Island Scrub Jays
}
\description{
Function to simulate open distance sampling data for the Island Scrub Jays, based on Sollmann et al (2015).

To recreate the data sets used in the book with R 3.6.0 or later, include \code{sample.kind="Rounding"} in the call to \code{set.seed}. This should only be used for reproduction of old results.
}
\usage{
issj.sim(B, db, lam, sigma, phi, gamma, npoints, nyrs, nbsize = -1.02)
}
\arguments{
  \item{B}{
Radius of the circle sampled; a site is a circle of radius B around a point.
}
  \item{db}{
Distance categories; a vector of cut points from 0 to B inclusive.
}
  \item{lam}{
Expected abundance per site, a vector of length \code{nsites} or an \code{nsites x 1} matrix.
}
  \item{sigma}{
Scale parameter of the half-normal detection function at each site,  a vector of length \code{nsites} or an \code{nsites x 1} matrix.
}
  \item{phi}{
Survival probability
}
  \item{gamma}{
Recruitment rate
}
  \item{npoints}{
Number of sites where point counts are conducted.
}
  \item{nyrs}{
Number of years
}
  \item{nbsize}{
Size parameter for the negative binomial distribution used to generate individual counts per site for year 1.
}
}
\value{
A list with the following elements:
 \item{NcList }{A list with one element per year, with distances of all animals from the point.}
 \item{detList }{A list with one element per year, a \code{npoints x nbands} matrix with the number of animals detected in each distance category for each point.}
 \item{N }{The (true) number of animals at each point for each year, a \code{nsites x nyears} matrix.}
 \item{cell }{The site IDs where point counts are conducted.}
 \item{y }{\code{npoints x nyears} matrix of total detections}
 \item{dclass }{a vector with the distance class for each animal detected}
 \item{site }{a corresponding vector with the site for each animal detected}
 \item{nsite }{the number of sites in the study area}
 \item{lam, phi, gamma, sigma }{the values of the corresponding arguments}
}
\references{
Sollmann, R., Gardner, B., Chandler, R.B., Royle, J.A., Sillett, T.S. (2015) An open population hierarchical distance sampling model. \emph{Ecology} 96, 325-331.

Kéry, M. & Royle, J.A. (2016) \emph{Applied Hierarchical Modeling in Ecology} AHM1 - 9.7.1.
}
\author{
Marc Kéry & Andy Royle, based on Sollmann et al (2015)
}

\examples{
# A toy example with just 20 sites
set.seed(2015)
tmp <- issj.sim(B = 300,
    db = c(0,50, 100, 150, 200, 250, 300),
    lam = c(3.01, 7.42, 20.51, 1.60, 0.42, 3.42, 8.24, 0.66, 0.32, 0.39, 0.46, 0.52,
      0.63, 0.36, 4.93, 0.47, 2.07, 0.42, 0.48, 0.47),
    sigma = c(110, 91, 70, 114, 135, 101, 88, 130, 133, 134, 134, 135, 131, 135, 100,
      135, 110, 135, 134, 135),
    phi = 0.6, gamma = 0.35,
    npoints = 15, nyrs = 4)
str(tmp)
# Compare the number detected with the true numbers present
with(tmp, cbind(y, N[cell, ]))
}
