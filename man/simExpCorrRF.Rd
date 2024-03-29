\name{simExpCorrRF}
\alias{simExpCorrRF}
\encoding{UTF-8}

\title{
Create a Gaussian random field with negative exponential correlation
}
\description{
Function creates Gaussian random field with negative exponential correlation and visualizes correlation and random field.

Previous versions used \pkg{RandomFields}, but that is not currently available on CRAN. If not available, \pkg{fields} will be used, but it cannot deal with large values of \code{size} and \code{theta}.
}
\usage{
simExpCorrRF(variance = 1, theta = 1, size = 50, seed = NA, show.plots = TRUE)
}

\arguments{
  \item{variance}{
variance of field.
}
  \item{theta}{
parameter governing spatial correlation (=1/phi) ("large theta means high correlation") Note that RMexp (which is used internally) is specified in terms of phi = 1/theta.
}
  \item{size}{
number of pixels along each side of the square site.
}
  \item{seed}{
used to generate reproducible output.
}
  \item{show.plots}{
if TRUE, the result will be displayed.
}
}

\value{
A list with the values of the input arguments and the following additional elements:
 \item{field }{the random field variable, a vector of length size^2}
 \item{grid }{the grid corresponding to field}
}
\references{
Kéry, M. & Royle, J.A. (2021) \emph{Applied Hierarchical Modeling in Ecology} AHM2 - 9.
}
\author{
Marc Kéry & Andy Royle
}

\examples{
# Run the function with default values and look at the output
str(tmp <- simExpCorrRF())
}
