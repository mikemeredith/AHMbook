\name{simHDS}
\alias{simHDS}
\encoding{UTF-8}

\title{
Simulate data under hierarchical distance sampling protocol (line or point)
}
\description{
Function simulates hierarchical distance sampling (HDS) data under either a line or a point transect protocol.
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
A list with the values of the arguments input, plus the following additional elements:
  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Kéry & Royle (2016) \emph{Applied Hierarachical Modeling in Ecology} Section 8.5.1
}
\author{
Kéry & Royle
}

\examples{
# Simulate a data set with the default arguments and look at the structure of the output:
tmp <- simHDS()
str(tmp)
}
