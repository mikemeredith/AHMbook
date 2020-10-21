\name{playRN}
\alias{playRN}
\encoding{UTF-8}
\title{
Function to play Royle-Nichols (RN) model
}
\description{
Function generates replicated count data under the binomial N-mixture model of Royle (2004), then 'degrades' the counts to detection/nondetection and fits the Royle-Nichols (RN) model (Royle & Nichols 2003) using \pkg{unmarked} and estimates site-specific abundance.
}
\usage{
playRN(M = 267, J = 3, mean.abundance = 1, mean.detection = 0.3,
    show.plots = TRUE, verbose = TRUE)
}
\arguments{
  \item{M}{
The number of sites.
}
  \item{J}{
The number of visits to each site.
}
  \item{mean.abundance}{
Expected abundance at each site.
}
  \item{mean.detection}{
Expected detection at each survey at each site.
}
  \item{show.plots}{
choose whether to show plots or not. Set to FALSE when using function in simulations.
}
  \item{verbose}{
if FALSE, output to the console will be suppressed.
}
}

\value{
A list with the following elements:
\item{nsites }{The number of sites, equal to \code{M} in the input}
\item{nvisits }{The number of visits, equal to \code{J} in the input}
\item{coef }{A named vector of coefficients for the linear models for expected number and detection probability}
\item{slope }{Slope of the regression of the estimated number on the true number; 1 if the model is perfect}
}
\references{
Royle, J.A. & Nichols, J.D. (2003) Estimating abundance from repeated presence-absence data or point counts, \emph{Ecology}, 84, 777-790.

Royle, J.A. (2004) N-mixture models for estimating population size from spatially replicated counts, \emph{Biometrics}, 60, 108-115.

Kéry, M. & Royle, J.A. (2016) \emph{Applied Hierarchical Modeling in Ecology} AHM1 - 6.13.1.
}
\author{
Marc Kéry & Andy Royle
}

\examples{
# Run a simulation with the default arguments and look at the results:
playRN()

\donttest{
# Execute the function using various settings
playRN(M = 100, J = 3, mean.abundance = 0.1)  # Increasing abundance
playRN(M = 100, J = 3, mean.abundance = 1)
playRN(M = 100, J = 3, mean.abundance = 5)
playRN(M = 100, J = 3, mean.detection = 0.3)  # Increasing detection
playRN(M = 100, J = 3, mean.detection = 0.5)
playRN(M = 100, J = 3, mean.detection = 0.7)
playRN(M = 100, J = 20)                       # More visits
playRN(M = 1000, J = 3)                       # More sites
}
}
