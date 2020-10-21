\name{simSpatialDSte}
\alias{simSpatialDSte}
\encoding{UTF-8}

\title{
Simulate data for replicate line transect surveys with temporary emigration
}
\description{
This simulates line transect distance sampling data with a spatial distribution of objects in a heterogeneous landscape where the parameter beta controls the effect of habitat. Multiple sample occasions are simulated and temporary emigration is allowed (parameter phi). Habitat is simulated according to a Gaussian random field model defined within the function. Uses a half normal detection model (if perp = TRUE) or a Gaussian hazard model (perp = FALSE).

To recreate the data sets used in the book with R 3.6.0 or later, include \code{sample.kind="Rounding"} in the call to \code{set.seed}. This should only be used for reproduction of old results.
}
\usage{
simSpatialDSte(nsites=28, dim=10, beta=1, lam0=2.5, nsurveys=4, sigma=3,
  phi=0.6, theta=2, show.plots=3)
}

\arguments{
  \item{nsites}{
number of sites
}
  \item{dim}{
number of pixels along each side of the square site
}
  \item{beta}{
the effect of habitat on the number of individuals in a pixel.
}
  \item{lam0}{
expected population size at each site
}
  \item{nsurveys}{
the number of replicate surveys
}
  \item{sigma}{
scale of half-normal detection function
}
  \item{phi}{
probability an individual is available for detection, ie, not temporarily emigrated.
}
  \item{theta}{
exponential correlation in the spatial covariate.
}
  \item{show.plots}{
the number of sites for which plots should be displayed; set to 0 to suppress plotting.
}
}

\value{
A list with the values of the input arguments and the following additional elements:
 \item{npixels }{the number of pixels in each site (= dim^2)}
 \item{B }{distance from line to edge of square (= dim/2)}
 \item{M }{true number of individuals at each site}
 \item{d }{perpendicular distance of each pixel from the line}
 \item{Habitat }{pixels x sites matrix, value of habitat covariate for each pixel}
 \item{y }{sites x pixels x surveys array, number of animals detected}
 \item{Counts }{sites x surveys matrix, number of animals detected (summed over pixels)}
}
\references{
Kéry, M. & Royle, J.A. (2016) \emph{Applied Hierarchical Modeling in Ecology}.

Mizel, J.D., Schmidt, J.H., & Lindberg, M.S. (2018) Accommodating temporary emigration in spatial distance sampling models. \emph{Journal of Applied Ecology}, 55, 1456-1464.

Kéry, M. & Royle, J.A. (2021) \emph{Applied Hierarchical Modeling in Ecology} AHM2 - 11.
}
\author{
Mizel et al (2018) Appendix S1, based in turn on Kéry & Royle (2016).
}

\examples{
# Run the function with default values and look at the output
str(tmp <- simSpatialDSte())
}
