\name{simpleNmix}
\alias{simpleNmix}
\encoding{UTF-8}
\title{
Simulate N-mixture data under a time-for-space substitution design
}
\description{
A simple function to simulate data under binomial N-mixture model where you have a single site that is surveyed over 'nyears' primary sampling periods ('seasons', 'years'), within each of which there are 'nreps' secondary samples.
}
\usage{
simpleNmix(nyears = 12, nreps = 4, beta0 = 2, beta1 = 0.1,
  alpha0 = 0.5, alpha1 = -0.1, alpha2 = 1, show.plot = TRUE)
}

\arguments{
  \item{nyears}{
Number of primary sampling periods.
}
  \item{nreps}{
Number of secondary samples within each primary period.
}
  \item{beta0}{
the intercept of a log-linear model of expected abundance (lambda).
}
  \item{beta1}{
the Time coefficient of a log-linear model for lambda.
}
  \item{alpha0}{
the intercept of a logit-linear model for detection (p).
}
  \item{alpha1}{
the Time coefficient of a logit-linear model for detection (p).
}
  \item{alpha2}{
the coefficient of a survey-specific covariate such as temperature (temp).
}
  \item{show.plot}{
choose whether to show plots or not. Set to FALSE when using function in simulations.
}
}

\value{
A list with the values of the arguments input and the following additional elements:
\item{N }{The realized number of individuals at each primary season, a vector of length \code{nyears}}
\item{C }{The number of individuals counted at each survey, a \code{nyears x nreps} matrix.}
\item{Time }{The Time covariate, a vector of length \code{nyears}.}
\item{temp }{The temperature covariate, a \code{nyears x nreps} matrix.}
\item{p }{The probability of detection, a \code{nyears x nreps} matrix.}
}
\references{
Kéry, M. & Royle, J.A. (2016) \emph{Applied Hierarchical Modeling in Ecology} AHM2 - 6.12.
}
\author{
Marc Kéry & Andy Royle
}

\examples{
# Simulate a data set with the default arguments and look at the structure of the output
tmp <- simpleNmix()
str(tmp)

}
