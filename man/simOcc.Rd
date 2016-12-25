\name{simOcc}
\alias{simOcc}
\encoding{UTF-8}
\title{
Simulate data for static occupancy models under wide range of conditions
}
\description{
Function to simulate occupancy measurements replicated at M sites during J occasions. Population closure is assumed for each site. Expected occurrence may be affected by elevation (elev), forest cover (forest) and their interaction. Expected detection probability may be affected by elevation, wind speed (wind) and their interaction.
}
\usage{
simOcc(M = 267, J = 3, mean.occupancy = 0.6, beta1 = -2, beta2 = 2, beta3 = 1,
  mean.detection = 0.3, time.effects = c(-1, 1),
  alpha1 = -1, alpha2 = -3, alpha3 = 0, sd.lp = 0.5,
  b = 2, show.plot = TRUE)
}
\arguments{
  \item{M}{
Number of spatial replicates (sites)
}
  \item{J}{
Number of temporal replicates (occasions)
}
  \item{mean.occupancy}{
Mean occurrence at value 0 of occurrence covariates
}
  \item{beta1}{
Main effect of elevation on occurrence
}
  \item{beta2}{
Main effect of forest cover on occurrence
}
  \item{beta3}{
Interaction effect on occurrence of elevation and forest cover
}
  \item{mean.detection}{
Mean detection prob. at value 0 of detection covariates
}
  \item{time.effects}{
bounds (on logit scale) for uniform distribution from which time effects gamma will be drawn
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
  \item{sd.lp}{
standard deviation of random site effects (on logit scale)
}
  \item{b}{
constant value of 'behavioural response' leading to 'trap-happiness' (if b > 0) or 'trap shyness' (if b < 0)
}
  \item{show.plot}{
if TRUE, plots of the data will be displayed; set to FALSE if you are running many simulations
}
}

\value{
A list with the values of the  input arguments and the following additional elements:

 \item{gamma }{The time effects, a vector of length J}
 \item{eps }{Individual random effects, a vector of length M}
 \item{elev }{Elevation, a vector of length M}
 \item{forest }{Forest cover, a vector of length M}
 \item{wind }{wind speed, a M x J matrix}
 \item{psi }{Probability of occurence, a vector of length M}
 \item{z }{Realised occurence (0/1), a vector of length M}
 \item{p }{probability of capture, possibly with a behavioural effect, a M x J matrix}
 \item{p0 }{probability of capture when not captured on previous occasion, a M x J matrix}
 \item{p1 }{probability of capture when captured on previous occasion, a M x J matrix}
 \item{y }{simulated capture history, a M x J matrix}
 \item{sumZ }{True number of occupied sites}
 \item{sumZ.obs }{Number of sites observed to be occupied}
 \item{psi.fs.true }{True proportion of occupied sites in sample (sumZ/N)}
 \item{psi.fs.obs }{Proportion of sites observed to be occupied (sumZ.obs/N)}
}
\references{
Kéry & Royle (2016) \emph{Applied Hierarachical Modeling in Ecology} Section 10.8

}
\author{
Marc Kéry
}

\examples{
# Generate data with the default arguments and look at the structure:
tmp <- simOcc()
str(tmp)
}
