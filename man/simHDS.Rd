\name{simHDS}
\alias{simHDS}
\encoding{UTF-8}

\title{
Simulate data under hierarchical distance sampling protocol (line or point)
}
\description{
The function simulates hierarchical distance sampling (HDS) data under either a line or a point transect protocol. At each site, it works with a strip of width \code{B*2} (for line transects) or a circle of radius \code{B} inscribed in a square of side \code{B*2} (for point transects).

The state process is simulated by first drawing a covariate value, "habitat", for each site from a Normal(0, 1) distribution. This is used in a log-linear regression with arguments \code{mean.lambda} and \code{beta.lam} to calculate the expected number of animals, lambda, in each strip or square.

For line transects, the distance from the line is drawn from a Uniform(0, B) distribution. For point transects, the animals are distributed randomly over the square before calculating the distance of each from the point. Observations of animals further than B from the point are discarded.

A detection covariate, "wind", for each site is drawn from a Uniform(-2, 2) distribution. This is used in a log-linear regression with arguments \code{mean.sigma} and \code{beta.sig} to calculate the scale parameter, sigma, of the half-normal detection function. Detections are simulated as Bernoulli trials with probability of success decreasing with distance from the line or point.
}
\usage{
simHDS(type=c("line", "point"), nsites = 100, mean.lambda = 2, beta.lam = 1,
  mean.sigma = 1, beta.sig = -0.5, B = 3, discard0 = TRUE, show.plot = TRUE)
}
\arguments{
  \item{type}{
type of transect, "line" or "point".
}
  \item{nsites}{
Number of sites (spatial replication)
}
  \item{mean.lambda}{
the expected value of lambda when the habitat covariate = 0; the intercept of the log-linear regression for lambda is log(mean.lambda).
}
  \item{beta.lam}{
the slope of the log-linear regression for lambda on a habitat covariate.
}
  \item{mean.sigma}{
the expected value of the scale parameter of the half-normal detection function when the wind speed = 0; the intercept of the log-linear regression for sigma is log(mean.sigma).
}
  \item{beta.sig}{
the slope of log-linear regression of scale parameter of the half-normal detection function on wind speed
}
  \item{B}{
the strip half-width or circle radius
}
  \item{discard0}{
If TRUE, subset to sites at which individuals were captured. You may or may not want to do this depending on how the model is formulated so be careful.
}
  \item{show.plot}{
choose whether to show plots or not. Set to FALSE when using function in simulations.
}
}

\value{
A list with the values of the arguments entered and the following additional elements:

\item{data }{simulated distance sampling data: a matrix with a row for each individual detected and 5 columns: site ID, status (1 if captured), x and y coordinates (NA for line transects), distance from the line or point; if \code{discard0 = FALSE}, sites with no detections will appear in the matrix with NAs in columns 2 to 5.}
\item{habitat }{simulated habitat covariate}
\item{wind }{simulated detection covariate}
\item{N }{simulated number of individuals at each site}
\item{N.true }{for point counts, the simulated number of individuals within the circle sampled}

}
\references{
Kéry, M. & Royle, J.A. (2016) \emph{Applied Hierarchical Modeling in Ecology} AHM1 - 8.5.1.
}
\author{
Marc Kéry & Andy Royle
}

\examples{
# Simulate a data set with the default arguments and look at the structure of the output
set.seed(123)
tmp <- simHDS()
str(tmp)
head(tmp$data)

tmp <- simHDS("point", discard0=FALSE)
str(tmp)
head(tmp$data, 10)
}
