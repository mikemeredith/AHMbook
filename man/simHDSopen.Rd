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
strip half width, or maximum distance
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
Wait...
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Kéry & Royle (2016) \emph{Applied Hierarachical Modeling in Ecology} Section 9.5.4.1
}
\author{
Marc Kéry & Andy Royle
}
\note{
For "point" the realized density is [(area of circle) /(area of  square)]*lambda
}

\examples{
tmp <- simHDSopen() # Generate data with default parameters
str(tmp)
}
