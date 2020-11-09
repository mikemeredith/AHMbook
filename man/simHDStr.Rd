\name{simHDStr}
\alias{simHDStr}
\encoding{UTF-8}
\title{
Simulate HDS time-removal or double-observer data
}
\description{
A general function for simulating hierarchical distance sampling (HDS) data combined with a time-removal (with 3 removal periods) or double-observer protocol, either for a line or a point transect protocol and with method = "removal" or method = "double". Also produces plots of the output.

At each site, it works with a strip of width \code{B*2} (for line transects) or a circle of radius \code{B} (for point transects).

The state process is simulated by first drawing a covariate value, "habitat", for each site from a Normal(0, 1) distribution. This is used in a log-linear regression with arguments \code{beta0} and \code{beta1} to calculate the expected number of groups in each strip or circle. Group size is simulated by first drawing from a Poisson distribution with parameter \code{lambda.group} then adding 1.

For line transects, the distance from the line is drawn from a Uniform(0, B) distribution. For point transects, the distance from the point is simulated from B*sqrt(Uniform(0,1)), which ensures a uniform distribution over the circle.

The group size is used in a log-linear regression with arguments \code{alpha0} and \code{alpha1} to calculate the scale parameter, sigma, of the half-normal detection function. Detections are simulated according to the selected protocol.

To recreate the data sets used in the book with R 3.6.0 or later, include \code{sample.kind="Rounding"} in the call to \code{set.seed}. This should only be used for reproduction of old results.
}
\usage{
simHDStr(type = c("line", "point"), method=c("removal", "double"), nsites = 200,
  lambda.group = 1, alpha0 = 0, alpha1 = 0, beta0 = 1, beta1 = 0.5,
  p.avail = 0.75, K = 3, p.double = c(0.4, 0.6),
  B = 3, discard0 = FALSE, show.plot = TRUE)
}
\arguments{
  \item{type}{
The type of distance transect, either "line" or "point".
}
  \item{method}{
Is the method time-removal ("removal") or double-observer ("double")?
}
  \item{nsites}{
Number of sites (spatial replication)
}
  \item{lambda.group}{
Poisson mean of group size
}
  \item{alpha0}{
intercept of log-linear model relating sigma of half-normal detection function to group size
}
  \item{alpha1}{
slope of log-linear model relating sigma of half-normal detection function to group size
}
  \item{beta0}{
intercept of log-linear model relating the Poisson mean of the number of groups per unit area to habitat
}
  \item{beta1}{
slope of log-linear model relating the Poisson mean of the number of groups per unit area to habitat
}
  \item{p.avail}{
overall availability probability (phi in text)
}
  \item{K}{
number of removal periods (of equal length)
}
  \item{p.double}{
detection probability for first and second observer
}
  \item{B}{
strip half width or radius of the circle
}
  \item{discard0}{
whether to discard or keep the data from sites with nobody detected
}
  \item{show.plot}{
choose whether to show plots or not. Set to FALSE when running simulations.
}
}

\value{
A list with the values of the arguments and the following additional elements:
\item{data }{simulated distance sampling data: a matrix with a row for each group detected and 7 columns: site ID, status (1 if captured), x and y coordinates (NA for line transects), distance from the line or point, group size, the time interval of removal or capture history; if \code{discard0 = FALSE}, sites with no detections will appear in the matrix with NAs in columns 2 to 7.}
\item{habitat }{simulated habitat covariate}
\item{M }{simulated number of groups at each site}
\item{M.true }{for point counts, the simulated number of groups within the circle sampled}
\item{params }{a vector with the input arguments}
}
\author{
Marc Kéry & Andy Royle
}
\references{
Kéry, M. & Royle, J.A. (2016) \emph{Applied Hierarchical Modeling in Ecology} AHM1 - 9.3.2 (time removal) and 9.4.1 (double observer).
}

\examples{
# Run with the default arguments and look at the structure of the output:
set.seed(123)
tmp <- simHDStr() # default: line, removal
str(tmp)
head(tmp$data)
tmp <- simHDStr("point", method="double")
str(tmp)
head(tmp$data)
}
