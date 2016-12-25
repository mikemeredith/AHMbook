\name{data.fn}
\alias{data.fn}
\encoding{UTF-8}
\title{
Simulate data for binomial mixture model
}
\description{
Function to simulate point counts replicated at M sites during J occasions. Population closure is assumed for each site. Expected abundance may be affected by elevation (elev), forest cover (forest) and their interaction. Expected detection probability may be affected by elevation, wind speed (wind) and their interaction.
}
\usage{
data.fn(M = 267, J = 3, mean.lambda = 2, beta1 = -2, beta2 = 2, beta3 = 1,
   mean.detection = 0.3, alpha1 = 1, alpha2 = -3, alpha3 = 0, show.plot = TRUE)
}
\arguments{
  \item{M}{
Number of spatial replicates (sites)
}
  \item{J}{
 Number of temporal replicates (occasions)
}
  \item{mean.lambda}{
Mean abundance at value 0 of abundance covariates
}
  \item{beta1}{
Main effect of elevation on abundance
}
  \item{beta2}{
Main effect of forest cover on abundance
}
  \item{beta3}{
Interaction effect on abundance of elevation and forest cover
}
  \item{mean.detection}{
Mean detection prob. at value 0 of detection covariates
}
  \item{alpha1}{
 Main effect of elevation on detection probability
}
  \item{alpha2}{
Main effect of wind speed on detection probability
}
  \item{alpha3}{
Interaction effect on detection of elevation and wind speed
}
  \item{show.plot}{
if TRUE, plots of the data will be displayed; set to FALSE if you are running many simulations
}
}

\value{
A list with the input arguments and the following additional elements:

 \item{elev }{Scaled elevation (a vector of length M)}
 \item{forest }{Scaled forest cover (a vector of length M)}
 \item{wind }{Scaled wind speed (an M x J matrix)}
 \item{lambda }{Expected number of individuals at each site (a vector of length M)}
 \item{N }{Realised number of individuals at each site (a vector of length M)}
 \item{p }{Probability of detection for each survey (an M x J matrix)}
 \item{C }{The number of detections for each survey (an M x J matrix)}
 \item{Ntotal }{Total abundance, \code{sum(N)} }
 \item{psi.true }{True occupancy in sample}
 \item{summaxC }{Sum of max counts (all sites)}
 \item{psi.obs }{ Observed occupancy in sample}

}
\note{
The colours used for points in some of the plots indicate different temporal replicates.
}
\references{
Kéry & Royle (2016) \emph{Applied Hierarachical Modeling in Ecology} Section 4.3
}
\author{
Marc Kéry and Andy Royle
}

\examples{
# Generate a simulated data set with default arguments and look at the structure:
tmp <- data.fn()
str(tmp)
}

