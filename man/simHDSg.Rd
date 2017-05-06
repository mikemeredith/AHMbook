\name{simHDSg}
\alias{simHDSg}
\encoding{UTF-8}

\title{
Simulate data under HDS protocol with groups
}
\description{
Simulates hierarchical distance sampling (HDS) data for groups under either a line or a point transect protocol and using a half-normal detection function (Buckland et al. 2001).
}
\usage{
simHDSg(type = c("line", "point"), nsites = 100, lambda.group = 0.75,
  alpha0 = 0, alpha1 = 0.5,
  beta0 = 1, beta1 = 0.5, B = 4, discard0 = TRUE, show.plot = TRUE)
}
\arguments{
  \item{type}{
The type of distance transect, either "line" or "point".
}
  \item{nsites}{
Number of sites (spatial replication)
}
  \item{lambda.group}{
Poisson mean of group size
}
  \item{alpha0}{
intercept of log-linear model relating sigma of the half-normal detection function to group size
}
  \item{alpha1}{
slope of log-linear model relating sigma of the half-normal detection function to group size
}
  \item{beta0}{
intercept of log-linear model relating the Poisson mean of the number of groups per unit area to habitat
}
  \item{beta1}{
slope of log-linear model relating the Poisson mean of the number of groups per unit area to habitat
}
  \item{B}{
strip half width or the radius of the circle
}
  \item{discard0}{
whether to discard or keep the data from sites with nobody detected
}
  \item{show.plot}{
choose whether to show plots or not. Set to FALSE when using function in simulations.
}
}
\value{
A list with the values of the arguments entered and the following additional elements:

\item{data }{simulated distance sampling data: a matrix with a row for each group detected and 6 columns: site ID, status (1 if captured), x and y coordinates (NA for line transects), distance from the line or point, group size; if \code{discard0 = FALSE}, sites with no detections will appear in the matrix with NAs in columns 2 to 6.}
\item{habitat }{simulated habitat covariate}
\item{N }{simulated number of groups at each site}
\item{N.true }{for point counts, the simulated number of groups within the circle sampled}
\item{groupsize }{group size for each of the groups observed}

}
\references{
Kéry & Royle (2016) \emph{Applied Hierarachical Modeling in Ecology} Section 9.2.1

Buckland et al (2001) \emph{Introduction to distance sampling: estimating abundance of biological populations}. Oxford University Press, Oxford, UK.
}

\examples{
# Run with the default arguments and look at the structure of the output:
set.seed(123)
tmp <- simHDSg()
str(tmp)
head(tmp$data)

str(simHDSg(type = "line"))     # Defaults for line transect data
str(simHDSg(type = "point"))    # Defaults for point transect data
str(simHDSg(lambda.group = 5))  # Much larger groups
str(simHDSg(lambda.group = 5, alpha1 = 0)) # No effect of groups size on p
}
