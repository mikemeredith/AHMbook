\encoding{UTF-8}
\name{simHDSopen}
\alias{simHDSopen}
\title{
Simulate open hierarchical distance sampling data
}
\description{
Simulates distance sampling data from a multi-season (or multi-year) model, incorporating habitat and detection covariates, temporary emigration, and a trend in abundance or density.
}
\usage{
simHDSopen(type=c("line", "point"), nsites = 100,
  mean.lam = 2, beta.lam = 0, mean.sig = 1, beta.sig = 0,
  B = 3, discard0 = TRUE, nreps = 2, phi = 0.7, nyears = 5, beta.trend = 0)
}
\arguments{
  \item{type}{
the transect protocol, either "line" or "point" .
}
  \item{nsites}{
Number of sites (spatial replication)
}
  \item{mean.lam}{
intercept of log-linear regression of expected lambda on a habitat covariate
}
  \item{beta.lam}{
slope of log-linear regression of expected lambda on a habitat covariate
}
  \item{mean.sig}{
intercept of log-linear regression of scale parameter of half-normal detection function on wind speed
}
  \item{beta.sig}{
slope of log-linear regression of scale parameter of half-normal detection function on wind speed
}
  \item{B}{
strip half width, or maximum distance from the observer for point counts
}
  \item{discard0}{
Discard sites at which no individuals were captured. You may or may not want to do this depending on how the model is formulated so be careful.
}
  \item{nreps}{
the number of distance sampling surveys within a period of closure in a season (or year)
}
  \item{phi}{
the availability parameter
}
  \item{nyears}{
the number of seasons (typically years)
}
  \item{beta.trend}{
loglinear trend of annual population size or density
}
}

\value{
A list with the values of the arguments entered and the following additional elements:

\item{data }{simulated distance sampling data: a list with a component for each year, each itself a list with a component for each replicate; this is a matrix with a row for each individual detected and 5 columns: site ID, status (1 if captured), x and y coordinates (NA for line transects), distance from the line or point; if \code{discard0 = FALSE}, sites with no detections will appear in the matrix with NAs in columns 2 to 5.}
\item{habitat }{simulated habitat covariate, a vector of length \code{nsites}}
\item{wind }{simulated detection covariate, a \code{nsites x nreps x nyears} array}
\item{M.true }{simulated number of individuals, a \code{nsites x nyears} matrix}
\item{K }{\code{= nreps}}
\item{Na }{the number of individuals available for detection, a \code{nsites x nreps x nyears} array}
\item{Na.real }{for point counts, the number of individuals available for detection within the circle sampled, a \code{nsites x nreps x nyears} array}

}
\references{
Kéry, M. & Royle, J.A. (2016) \emph{Applied Hierarchical Modeling in Ecology} Section 9.5.4.1
}
\author{
Marc Kéry & Andy Royle
}
\note{
For "point" the realized density is [(area of circle) /(area of  square)]*lambda
}

\examples{
set.seed(123)
tmp <- simHDSopen() # Generate data with default parameters
str(tmp)
head(tmp$data[[1]][[1]])

tmp <- simHDSopen("point")
str(tmp)
head(tmp$data[[1]][[1]])

tmp <- simHDSopen(discard0=FALSE)
str(tmp)
head(tmp$data[[1]][[1]])
}
