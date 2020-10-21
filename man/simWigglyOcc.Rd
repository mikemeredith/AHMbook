\name{wigglyOcc}
\alias{wigglyOcc}
\encoding{UTF-8}
\title{
Simulate static occupancy data
}
\description{
Function to generate a static occupancy data set with really wiggly covariate relationships
in occupancy and detection probability.

To recreate the data sets used in the book with R 3.6.0 or later, include \code{sample.kind="Rounding"} in the call to \code{set.seed}. This should only be used for reproduction of old results.
}
\usage{
wigglyOcc(seed = 1, show.plot = TRUE, verbose = TRUE)
}
\arguments{
  \item{seed}{
Seed for random number generator
}
  \item{show.plot}{
choose whether to show plots or not. Set to FALSE when using function in simulations.
}
  \item{verbose}{
if TRUE, output will be written to the console.
}
}

\value{
A list with the following elements:
  \item{M }{Number of sites}
  \item{J }{Number of replicate surveys}
  \item{Xsite }{Simulated site covariate, a vector of length M}
  \item{Xsurvey }{Simulated survey covariate, a M x J matrix}
  \item{psi }{Occupancy probability, a vector of length M}
  \item{z }{Realized occupancy, a 0/1 vector of length M}
  \item{p }{Detection probability, a M x J matrix}
  \item{y }{detection history, a M x J matrix of 0/1}
}
\references{
Kéry, M. & Royle, J.A. (2016) \emph{Applied Hierarchical Modeling in Ecology} AHM1 - 10.14.

}
\author{
Marc Kéry & Andy Royle
}

\examples{
# Generate data with the default arguments and look at the structure:
tmp <- wigglyOcc()
str(tmp)
}
