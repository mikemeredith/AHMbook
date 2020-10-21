\name{simMultMix}
\alias{simMultMix}
\encoding{UTF-8}

\title{
Simulate the multinomial-mixture model
}
\description{
Simulation of "removal" data from a multinomial-mixture model.
}
\usage{
simMultMix(nsites = 100, nsurveys = 3, nyears = 4,
  lambda = 3, theta = 0.5, p = 0.3)
}

\arguments{
  \item{nsites}{
number of sites.
}
  \item{nsurveys}{
number of replicate (secondary) samples
}
  \item{nyears}{
number of primary samples: years, seasons etc.
}
  \item{lambda}{
expected local population size.
}
  \item{theta}{
availability, the proportion of the population available for detection.
}
  \item{p}{
detection probability.
}
}

\value{
A list with the values of the input arguments and the following additional elements:
 \item{M }{true local population size, vector length nsites}
 \item{N }{true number available for detection, vector length nsites}
 \item{y }{number detected, nsites x nyears x nsurveys}
 \item{y2d }{number detected as a matrix, nsites x (nyears*nsurveys)}

}
\references{
Chandler, R.B., Royle, J.A. & King. D.I. (2011) Inference about density and temporary emigration in unmarked populations. \emph{Ecology}, 92, 1429-1435.

Kéry, M. & Royle, J.A. (2021) \emph{Applied Hierarchical Modeling in Ecology} AHM2 - 2.7.1.
}
\author{
Marc Kéry & Andy Royle, adapting code in documentation for \code{\link[unmarked]{gmultmix}} by Richard Chandler.
}

\examples{
# Simulate a data set with the default arguments and look at the structure of the output:
tmp <- simMultMix()
str(tmp)
}
