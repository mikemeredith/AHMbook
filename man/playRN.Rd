\name{playRN}
\alias{playRN}
\encoding{UTF-8}
\title{
Function to play Royle-Nichols model
}
\description{
Function generates replicated count data under the Nmix model of Royle (2004), then 'degrades' the data to detection/nondetection and fits the RN model (Royle & Nichols 2003) using \pkg{unmarked} and estimates site-specific abundance.
}
\usage{
playRN(M = 267, J = 3, mean.abundance = 1, mean.detection = 0.3, show.plot = TRUE)
}
\arguments{
  \item{M}{
The number of sites.
}
  \item{J}{
The number of visits to each site.
}
  \item{mean.abundance}{
Expected abundance at each site.
}
  \item{mean.detection}{
Expected detection at each survey at each site.
}
  \item{show.plot}{
choose whether to show plots or not. Set to FALSE when using function in simulations.
}
}

\value{
A list with the following elements:
\item{nsite }{The number of sites, equal to \code{M} in the input}
\item{nvisit }{The number of visits, equal to \code{J} in the input}
\item{coef }{A named vector of coefficients for the linear models for expected number and detection probability}
\item{slope }{Slope of the regression of the estimated number on the true number; 1 if the model is perfect}
}
\references{
Kéry & Royle (2016) \emph{Applied Hierarachical Modeling in Ecology} Section 6.13.1
}
\author{
Kéry & Royle
}

\examples{
# Run a simulation with the default arguments and look at the results:
playRN()
}

