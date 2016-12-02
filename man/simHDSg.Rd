\name{simHDSg}
\alias{simHDSg}
\encoding{UTF-8}

\title{
Simulate data under HDS protocol with groups
}
\description{
Simulates hierarchical distance sampling (HDS) data for groups under either a line (type = "line") or a point (type = "point") transect protocol and using a half-normal detection function (Buckland et al. 2001).
}
\usage{
simHDSg(type = "line", nsites = 100, lambda.group = 0.75, alpha0 = 0,
  alpha1 = 0.5, beta0 = 1, beta1 = 0.5, B = 4, discard0 = TRUE)
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
strip half width
}
  \item{discard0}{
whether to discard or keep the data from sites with nobody detected
}
}
\value{
A list with the values of the arguments entered and the following additional elements:

\item{data }{simulated distance sampling data}
\item{habitat }{simulated habitat covariate}
\item{N }{simulated number of groups at each site}
\item{N.true }{for point counts, the simulated number of groups within the circle sampled}
\item{groupsize }{group size for each of the oups observed}

}
\references{
Kéry & Royle (2016) \emph{Applied Hierarachical Modeling in Ecology} Section 9.2.1

Buckland et al (2001) 
}

\examples{
# Run with the default arguments and look at the structure of the output:
str(simHDSg)
}
