\name{sim3Occ}
\alias{sim3Occ}
\encoding{UTF-8}
\title{
Simulate data for static 3-level occupancy models
}
\description{
Function generates 3-level occupancy data with possibility of site-specific random variation at every level, "time effects" at the middle and the lower levels and effects of one distinct covariate at each level.
}
\usage{
sim3Occ(nunit = 100, nsubunit = 5, nrep = 3, mean.psi = 0.8, beta.Xpsi = 1, sd.logit.psi = 0, mean.theta = 0.6, theta.time.range = c(-1, 1), beta.Xtheta = 1, sd.logit.theta = 0, mean.p = 0.4, p.time.range = c(-2, 2), beta.Xp = -1, sd.logit.p = 0)
}
\arguments{
  \item{nunit}{
Number of main units (large quadrats)
}
  \item{nsubunit}{
Number of subunits (nested subsamples within each main unit)
}
  \item{nrep}{
Number of replicate surveys in every subunit
}
  \item{mean.psi}{
Mean large-scale, unit-level occupancy probability (psi)
}
  \item{beta.Xpsi}{
effect on psi of covariate A (at main unit level)
}
  \item{sd.logit.psi}{
SD of logit(psi), unstructured site variation in psi
}
  \item{mean.theta}{
Mean small-scale (subunit) occupancy probability (theta)
}
  \item{theta.time.range}{
range of theta 'intercepts' for subunits
}
  \item{beta.Xtheta}{
effect on theta of covariate B (at subunit level)
}
  \item{sd.logit.theta}{
SD of logit(theta), unstructured site variation in theta
}
  \item{mean.p}{
Mean per-survey detection probability
}
  \item{p.time.range}{
range of p 'intercepts' for replicates
}
  \item{beta.Xp}{
effect on p of covariate C (unit by subunit by replicate)
}
  \item{sd.logit.p}{
SD of logit(p)
}
}
\value{
A list with the values of the input arguments and the following additional elements:
 \item{theta.time.effect }{Simulated time effect on theta, a vector of length \code{nsubunit}}
 \item{p.time.effect }{Simulated time effect on p, a vector of length \code{nrep}}
 \item{p }{Detection probabiliy, a \code{nunit x nsubunit x nrep} array}
 \item{z }{Occupancy indicator for main units, a \code{nunit x 1} array of 0/1}
 \item{a }{Occupancy indicator for subunits, a \code{nunit x nsubunit} array of 0/1}
 \item{y }{Detection array, a \code{nunit x nsubunit x nrep} array of 0/1}
 \item{sum.z }{True number of occupied main units}
 \item{obs.sum.z }{Observed number of occupied main units}
 \item{sum.z.a }{Number of units with >=1 occupied, surveyed subunit}
 \item{covA }{Simulated covariate A, a vector of length \code{nunit}}
 \item{covB }{Simulated covariate B, a \code{nunit x nsubunit} matrix}
 \item{covC }{Simulated covariate C, a \code{nunit x nsubunit x nrep} array}
}
\references{
Kéry & Royle (2016) \emph{Applied Hierarachical Modeling in Ecology} Section 10.10

}
\author{
Marc Kéry
}

\examples{
# Generate data with the default arguments and look at the structure:
tmp <- sim3Occ()
str(tmp)
}
