\name{wigglyOcc}
\alias{wigglyOcc}
\encoding{UTF-8}
\title{
Simulate static occupancy data
}
\description{
Function to generate a static occupancy data set with really wiggly covariate relationships
in occupancy and detection probability}
\usage{
wigglyOcc(seed = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seed}{
Seed for random number generator
}
}

\value{
A list with the following elements:
  \item{M }{Number of sites}
  \item{J }{Number of replicate surveys}
  \item{Xsite }{Simulated site covariate, a vector of length M}
  \item{Xsurvey }{Simulated survey covariate, a M x J matrix}
  \item{x.index }{??, a vector of length M x J}
  \item{psi }{Occupancy probability, a vector of length M}
  \item{z }{Realised occupancy, a 0/1 vector of length M}
  \item{p }{Detection probability, a M x J matrix}
  \item{p.ordered }{??, a vector of length M x J}
  \item{y }{detection history, a M x J matrix of 0/1}
}
\references{
Kéry & Royle (2016) \emph{Applied Hierarachical Modeling in Ecology} Section 10.14

}
\author{
Kéry & Royle
}

\examples{
# Generate data with the default arguments and look at the structure:
tmp <- wigglyOcc()
str(tmp)
}
