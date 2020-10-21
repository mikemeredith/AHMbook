\name{sim3Occ}
\alias{sim3Occ}
\encoding{UTF-8}
\title{
Simulate detection/nondetection data for static 3-level occupancy models
}
\description{
Function generates 3-level occupancy data with possibility of site-specific random variation at every level, "time effects" at the middle and the lower levels and effects of one distinct covariate at each level.
}
\usage{
sim3Occ(nunits = 100, nsubunits = 5, nreps = 3,
  mean.psi = 0.8, beta.Xpsi = 1, sd.logit.psi = 0,
  mean.theta = 0.6, theta.time.range = c(-1, 1), beta.Xtheta = 1, sd.logit.theta = 0,
  mean.p = 0.4, p.time.range = c(-2, 2), beta.Xp = -1, sd.logit.p = 0,
  show.plot = TRUE, verbose = TRUE)
}
\arguments{
  \item{nunits}{
Number of main units (large quadrats)
}
  \item{nsubunits}{
Number of subunits (nested subsamples within each main unit)
}
  \item{nreps}{
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
  \item{show.plot}{
choose whether to show plots or not. Set to FALSE when using function in simulations.
}
  \item{verbose}{
if TRUE, output will be written to the console.
}
}
\value{
A list with the values of the input arguments and the following additional elements:
 \item{theta.time.effect }{Simulated time effect on theta, a vector of length \code{nsubunits}}
 \item{p.time.effect }{Simulated time effect on p, a vector of length \code{nreps}}
 \item{p }{Detection probabiliy, a \code{nunits x nsubunits x nreps} array}
 \item{z }{Occupancy indicator for main units, a \code{nunits x 1} array of 0/1}
 \item{a }{Occupancy indicator for subunits, a \code{nunits x nsubunits} array of 0/1}
 \item{y }{Detection array, a \code{nunits x nsubunits x nreps} array of 0/1}
 \item{sum.z }{True number of occupied main units}
 \item{obs.sum.z }{Observed number of occupied main units}
 \item{sum.z.a }{Number of units with >=1 occupied, surveyed subunit}
 \item{covA }{Simulated covariate A, a vector of length \code{nunits}}
 \item{covB }{Simulated covariate B, a \code{nunits x nsubunits} matrix}
 \item{covC }{Simulated covariate C, a \code{nunits x nsubunits x nreps} array}
}
\references{
Kéry, M. & Royle, J.A. (2016) \emph{Applied Hierarchical Modeling in Ecology} AHM1 - 10.10.
}

\author{
Marc Kéry & Andy Royle
}

\examples{
# Generate data with the default arguments and look at the structure:
tmp <- sim3Occ()
str(tmp)

# 'Null' model (model 1)
str(data <- sim3Occ(nunits = 100, nsubunits = 5, nreps = 3, mean.psi = 0.8,
  beta.Xpsi = 0, sd.logit.psi = 0, mean.theta = 0.6, theta.time.range = c(0, 0),
  beta.Xtheta = 0, sd.logit.theta = 0, mean.p = 0.4, p.time.range = c(0,0),
  beta.Xp = 0, sd.logit.p = 0))

# No covariate effects, no random variability (model 2)
str(data <- sim3Occ(nunits = 100, nsubunits = 5, nreps = 3, mean.psi = 0.8,
  beta.Xpsi = 0, sd.logit.psi = 0, mean.theta = 0.6, theta.time.range = c(-1, 1),
  beta.Xtheta = 0, sd.logit.theta = 0, mean.p = 0.4, p.time.range = c(-2,2),
  beta.Xp = 0, sd.logit.p = 0))

# All covariate effects, but no random variability (model 3)
str(data <- sim3Occ(nunits = 100, nsubunits = 5, nreps = 3, mean.psi = 0.8,
  beta.Xpsi = 1, sd.logit.psi = 0, mean.theta = 0.6, theta.time.range = c(-1, 1),
  beta.Xtheta = 1, sd.logit.theta = 0, mean.p = 0.4, p.time.range = c(-2,2),
  beta.Xp = -1, sd.logit.p = 0))

# Most complex model with all effects allowed for by sim function (model 4)
str(data <- sim3Occ(nunits = 100, nsubunits = 5, nreps = 3, mean.psi = 0.8,
  beta.Xpsi = 1, sd.logit.psi = 1, mean.theta = 0.6, theta.time.range = c(-1, 1),
  beta.Xtheta = 1, sd.logit.theta = 1, mean.p = 0.4, p.time.range = c(-2,2),
  beta.Xp = -1, sd.logit.p = 1))
}
