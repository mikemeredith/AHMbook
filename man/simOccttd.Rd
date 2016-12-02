\name{simOccttd}
\alias{simOccttd}
\encoding{UTF-8}
\title{
Simulate time-to-detection occupancy data
}
\description{
Function simulates time-to-detection occupancy design data under model of Garrard et al. (Austral Ecology, 2008), also see Bornand et al. (MEE, 2014)
}
\usage{
simOccttd(M = 250, mean.psi = 0.4, mean.lambda = 0.3, beta1 = 1, alpha1 = -1, Tmax = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{M}{
Number of sites
}
  \item{mean.psi}{
intercept of occupancy probability
}
  \item{mean.lambda}{
intercept of Poisson rate parameter
}
  \item{beta1}{
slope of continuous covariate B on logit(psi)
}
  \item{alpha1}{
slope of continuous covariate A on log(lambda)
}
  \item{Tmax}{
maximum search time (in arbitrary units, which are same as response), response will be censored at Tmax
}
}

\value{
A list with the values of the arguments input and the following additional elements:
 \item{covA }{Simulated values of covariate A, a vector of length M}
 \item{covB }{Simulated values of covariate B, a vector of length M}
 \item{psi }{Probability of occurence at each site, a vector of length M}
 \item{z }{Realised occurence at each site, a 0/1 vector of length M}
 \item{ttd.temp }{Uncensored simulated time-to-detection at each site, a vector of length M}
 \item{ttd }{Censored simulated time-to-detection at each site, a vector of length M}
 \item{d }{Censoring indicator, a 0/1 vector of length M}
 \item{sum.z }{Total number of sites occupied}
 \item{n.obs }{Total number of sites where the species was observed}

}
\references{
Kéry & Royle (2016) \emph{Applied Hierarachical Modeling in Ecology} Section 10.12.1

Garrard, G.E., Bekessy, S.A., McCarthy, M.A., & Wintle, B.A. (2008) When have we looked hard enough? A novel method for setting minimum survey effort protocols for flora surveys. \emph{Austral Ecology}, 33, 986-998.

Bornand, C.N., Kéry, M., Bueche, L., & Fischer, M. (2014) Hide and seek in vegetation: time-to-detection is an efficient design for estimating detectability and occurrence. \emph{Methods in Ecology and Evolution}, 5, 433-442.
}
\author{
Marc Kéry
}

\examples{
# Generate data with the default arguments and look at the structure:
tmp <- simOccttd()
str(tmp)
}
