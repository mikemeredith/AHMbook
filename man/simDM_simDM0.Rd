\name{simDM}
\alias{simDM}
\alias{simDM0}
\encoding{UTF-8}

\title{
Simulate count data under the dynamic N-mixture model of Dail-Madsen
}
\description{
Simulation for multiple-visit data (from \code{\link[unmarked]{pcountOpen}} help file in package \pkg{unmarked}). \code{simDM0} has no covariates while \code{simDM} allows for covariates. Both functions assume constant time intervals between primary periods.
}
\usage{
simDM0(nsites = 50, nsurveys = 3, nyears = 5,
  lambda = 4, gamma = 1.5, phi = 0.8, p = 0.7, show.plots=TRUE)

simDM(nsites = 50, nsurveys = 3, nyears = 5,
  mean.lambda = 4, mean.gamma.rel = 0.5,
  mean.phi = 0.8, mean.p = 0.7,
  beta.lam = 1, beta.gamma = 1, beta.phi = -1, beta.p = -1,
  show.plots=TRUE)
}

\arguments{
  \item{nsites}{
number of sites.
}
  \item{nsurveys}{
number of replicate (secondary) samples within period of closure.
}
  \item{nyears}{
number of primary samples: years, seasons etc.
}
  \item{lambda}{
Initial expected abundance.
}
  \item{gamma}{
recruitment rate.
}
  \item{phi}{
apparent survival rate.
}
  \item{p}{
detection probability.
}
  \item{mean.lambda}{
Initial expected abundance at cov.lam = 0.
}
  \item{mean.gamma.rel}{
recruitment rate at cov.gamma = 0.
}
  \item{mean.phi}{
apparent survival rate at cov.phi = 0.
}
  \item{mean.p}{
detection probability at cov.p = 0.
}
  \item{beta.lam}{
the slope of parameter lambda (link transformed) on the cov.lam covariate
}
  \item{beta.gamma}{
the slope of parameter gamma (link transformed) on the cov.gamma covariate
}
  \item{beta.phi}{
the slope of parameter phi (link transformed) on the cov.phi covariate
}
  \item{beta.p}{
the slope of parameter p (link transformed) on the cov.p covariate
}
  \item{show.plots}{
if TRUE, plots are produced. Set this to FALSE when running simulations.
}
}

\value{
For \code{simDM0}, a list with the values of the input arguments and the following additional elements:
 \item{N }{true number of individuals, nsites x nyears}
 \item{S }{number of survivors, nsites x (nyears-1)}
 \item{R }{number of recruits, nsites x (nyears-1)}
 \item{y }{number detected, nsites x nyears x nsurveys}
 \item{yy }{number detected as a matrix, nsites x (nyears*nsurveys)}

\code{simDM} has the following additional elements:
 \item{cov.lam }{covariate for lambda generated from Uniform(-1, 1), nsites vector}
 \item{cov.gamma }{covariate for gamma generated from Uniform(-1, 1), nsites vector}
 \item{cov.phi }{covariate for phi generated from Uniform(-1, 1), nsites vector}
 \item{cov.p }{covariate for p generated from Uniform(-1, 1), nsites x nyears x nsurveys}
 \item{ccov.p }{cov.p formatted as a matrix, nsites x (nyears*nsurveys)}
}
\references{
Dail, D. & Madsen, L. (2011) Models for estimating abundance from repeated counts of an open metapopulation. \emph{Biometrics}, 67, 577-587.

Kéry, M. & Royle, J.A. (2021) \emph{Applied Hierarchical Modeling in Ecology} AHM2 - 2.5.
}
\author{
Marc Kéry & Andy Royle, based in part on code written by Richard Chandler.
}

\examples{
# Simulate a data set with the default arguments and look at the structure of the output:
tmp0 <- simDM0()
str(tmp0)

tmp <- simDM()
str(tmp)

}
