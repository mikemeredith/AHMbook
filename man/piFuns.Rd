\name{piFuns}
\alias{piFuns}
\alias{instRemPiFun}
\alias{crPiFun}
\alias{crPiFun.Mb}
\alias{MhPiFun}
\encoding{UTF-8}
\title{
Compute multinomial cell probabilities
}
\description{
These functions are used internally in the multinomial-Poisson models \code{multinomPois} and \code{gmultmix} in package \pkg{unmarked} to calculate the multinomial cell probabilities from the "raw" detection probabilities. See \link[unmarked]{piFuns}. These functions are customized for the examples in chapter 7 of AHM1 and only allow for 3 capture occasions.

See \link[unmarked]{makePiFuns} for functions to construct piFuns with more occasions.
}
\usage{
instRemPiFun(p)

crPiFun(p)

crPiFun.Mb(p)

MhPiFun(p)
}
\arguments{
  \item{p}{
a numeric matrix with rows for each site and 3 columns, see Details.
}
}
\details{
\code{instRemPiFun} defines the relationship between the multinomial cell probabilities and the underlying detection probability parameters in a removal design with 3 sampling periods of unequal length, specifically 2, 3 and 5 minutes. The columns of \code{p} give the detection probabilities per unit time (ie, per minute) for each site and sampling period. See AHM1 Section 7.7. This is the same as the function defined in the Examples section of \link[unmarked]{piFuns} in \pkg{unmarked}.

\code{crPiFun} defines a pi function for capture-recapture design with 3 capture occasions. The columns of \code{p} give the detection probabilities for the three occasions. See AHM1 Section 7.8.5. \bold{NOTE} that this is not the same as the custom \code{crPiFun} defined in Section 7.9.1.

\code{crPiFun.Mb} defines a pi function for capture-recapture design with 3 capture occasions with a behavioral response. Detection probabilities do not vary with time. Column #1 of \code{p} gives the probability of detection for animals not previously detected; column #3 gives detection probability after the first detection; column #2 is ignored. See AHM1 Section 7.8.2.

\code{MhPiFun} defines a pi function for a model with individual detection heterogeneity, modeled as a random effect with a logit-normal distribution. Column #1 of \code{p} gives the mean probability of detection; \code{p[1, 2]} is a value in [0, 1] which controls the scale parameter for the normal distribution; other entries are ignored. See Section 7.8.3.
}
\value{
A matrix with a row for each site corresponding to the rows of \code{p} and...

...for \code{instRemPiFun}, a column for each detection occasion with the multinomial probability.

...for other functions, a column for each detection history; for 3 detection occasions these are 001, 010, 011, 100, 101, 110, 111.
}
\references{
Kéry, M. & Royle, J.A. (2016) \emph{Applied Hierarchical Modeling in Ecology} AHM1 - 7.
}
\author{
Richard Chandler, Andy Royle
}

\examples{
# Capture probs for 5 sites, with 3 survey periods
( p <- matrix(0.4, nrow=5, ncol=3) )
# Removal model
instRemPiFun(p)
# The corresponding obsToY matrix:
matrix(1, 2, 3)

# Capture-recapture model
crPiFun(p)
# The corresponding obsToY matrix:
matrix(1, 3, 7)

# Capture-recapture model with behavioural response
( pMb <- cbind(rep(0.7, 5), NA, 0.3) )
crPiFun.Mb(pMb)
# The corresponding obsToY matrix:
matrix(1, 3, 7)

# Capture-recapture model with heterogeneity
pMh <- cbind(rep(0.4, 5), NA, NA)
pMh[1, 2] <- 0.3
pMh
MhPiFun(pMh)
# The corresponding obsToY matrix:
matrix(1, 3, 7)
}
