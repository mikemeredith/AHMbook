\name{simpleNmix}
\alias{simpleNmix}
\encoding{UTF-8}
\title{
Simulate N-mixture data under a time-for-space substitution design
}
\description{
A simple function to simulate data under binomial N-mixture model where you have a single site that is survyed over 'nyear' primary sampling periods ('seasons', 'years'), within each of which there are 'nrep' secondary samples.
}
\usage{
simpleNmix(nyear = 12, nrep = 4, beta0 = 2, beta1 = 0.1, alpha0 = 0.5, alpha1 = -0.1, alpha2 = 1)
}

\arguments{
  \item{nyear}{
Number of primary sampling periods.
}
  \item{nrep}{
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
}

\value{
A list with the values of the arguments input and the following additional elements"
\item{N }{The realised number of individuals at each primary season, a vector of length \code{nyear}}
\item{C }{The number of individuals counted at each survey, a \code{nyear x nrep} matrix.}
\item{Time }{The Time covariate, a vector of length \code{nyear}.}
\item{temp }{The temperature covariate, a \code{nyear x nrep} matrix.}
\item{p }{The probability of detection, a \code{nyear x nrep} matrix.}
}
\references{
Kéry & Royle (2016) \emph{Applied Hierarachical Modeling in Ecology} Section 6.12
}
\author{
Kéry & Royle
}

\examples{
# Simulate a data set with the default arguments and look at the structure of the output:
tmp <- simpleNmix()
str(tmp)

}
