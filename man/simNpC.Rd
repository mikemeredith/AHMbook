\name{simNpC}
\alias{simNpC}
\encoding{UTF-8}

\title{
Generate counts from a single population observed over T years under a binomial observation process
}
\description{
Generates counts from a single population observed over T years and which can be observed with or without imperfect detection. The goal is to focus on what happens with relative-abundance inference when temporal patterns in abundance are confounded with temporal patterns in detection probability. Hence, we can simulate a stable population or one with linear increase or decrease with specified start and end points, and around which there is Poisson noise. The observed counts are Binomial outcomes with a detection probability which can similarly be chosen to be constant or change linearly over time.
}
\usage{
simNpC(T = 20, expN = c(100, 75), dp = c(0.5, 0.5), show.plot = TRUE)
}
\arguments{
  \item{T}{
The length of the time series.
}
  \item{expN}{
The expected abundance at start and end of period, linear trend.
}
  \item{dp}{
The detection probability at start and end of period, linear trend.
}

  \item{show.plot}{
Choose whether to show plots or not. Set to FALSE when using function in simulations.
}
}
\value{
A list with the values of the arguments entered and the following additional elements:

\item{lambda }{T vector, expected abundance for each year.}
\item{p }{T vector, detection probability (dp) for each year.}
\item{N }{T vector, realized abundance.}
\item{C }{T vector, observed counts.}
}
\references{
Kéry, M. & Royle, J.A. (2021) \emph{Applied Hierarchical Modeling in Ecology} AHM2 - 1.2.

}
\author{
Marc Kéry & Andy Royle
}

\examples{
# Run with the default arguments and look at the structure of the output:
set.seed(123)
tmp <- simNpC()
str(tmp)
tmp$C


}
