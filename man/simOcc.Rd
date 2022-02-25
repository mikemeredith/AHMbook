\name{simOcc}
\alias{simOcc}
\alias{simOccCat}
\encoding{UTF-8}
\title{
Simulate data for static occupancy models under wide range of conditions
}
\description{
Functions to simulate static occupancy measurements replicated at M sites during J occasions. Population closure is assumed for each site. Expected occurrence may be affected by elevation (elev), forest cover (forest) and their interaction. Expected detection probability may be affected by elevation, wind speed (wind) and their interaction.

\code{simOccCat} allows for a categorical habitat-type covariate (HAB) to affect occurrence and another, the observer (OBS), to affect detection.
}
\usage{
simOcc(M = 267, J = 3, mean.occupancy = 0.6, beta1 = -2, beta2 = 2, beta3 = 1,
  mean.detection = 0.3, time.effects = c(-1, 1),
  alpha1 = -1, alpha2 = -3, alpha3 = 0, sd.lp = 0.5,
  b = 2, show.plots = TRUE)

simOccCat(M = 267, J = 3, mean.occupancy = 0.6, beta1 = 0, beta2 = 0, beta3 = 0,
  mean.detection = 0.3, time.effects = c(0, 0),
  alpha1 = 0, alpha2 = 0, alpha3 = 0, sd.lp = 0, b = 0,
  nHab = 5, range.HAB = 2, nObs = 10, range.OBS = 4,
  show.plots = TRUE)

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
constant value of 'behavioral response' leading to 'trap-happiness' (if b > 0) or 'trap shyness' (if b < 0)
}
  \item{nHab}{
number of categories to simulate for the site-specific categorical covariate, HAB.
}
  \item{range.HAB}{
a non-negative number that controls the size of the HAB effect.
}
  \item{nObs}{
number of observers, the categories to simulate for the survey-specific categorical covariate, OBS.
}
  \item{range.OBS}{
a non-negative number that controls the size of the OBS effect.
}
  \item{show.plots}{
if TRUE, plots of the data will be displayed; set to FALSE if you are running simulations
}
}

\value{
A list with the values of the  input arguments and the following additional elements:

 \item{gamma }{The time effects, a vector of length J}
 \item{eps }{Individual random effects, a vector of length M}
 \item{elev }{Elevation, a vector of length M}
 \item{forest }{Forest cover, a vector of length M}
 \item{wind }{wind speed, a M x J matrix}
 \item{psi }{Probability of occurrence, a vector of length M}
 \item{z }{Realized occurrence (0/1), a vector of length M}
 \item{p }{probability of capture, possibly with a behavioral effect, a M x J matrix}
 \item{p0 }{probability of capture when not captured on previous occasion, a M x J matrix}
 \item{p1 }{probability of capture when captured on previous occasion, a M x J matrix}
 \item{y }{simulated capture history, a M x J matrix}
 \item{sumZ }{True number of occupied sites}
 \item{sumZ.obs }{Number of sites observed to be occupied}
 \item{psi.fs.true }{True proportion of occupied sites in sample (sumZ/N)}
 \item{psi.fs.obs }{Proportion of sites observed to be occupied (sumZ.obs/N)}

In addition, \code{simOccCat} has
 \item{HAB }{Categorical site covariate, a vector of length M}
 \item{OBS }{Categorical detection covariate, a M x J matrix}
 \item{coefHAB }{The effect of each category of HAB on occurrence}
 \item{coefOBS }{The effect of each category of OBS on detection}

}
\references{
Kéry, M. & Royle, J.A. (2016) \emph{Applied Hierarchical Modeling in Ecology} AHM1 - 10.8.

}
\author{
Marc Kéry, Andy Royle, Gesa von Hirschheydt, Mike Meredith,
}

\examples{
# Generate data with the default arguments and look at the structure
tmp <- simOcc()
str(tmp)

# Simplest possible occupancy model, with constant occupancy and detection
str(simOcc(mean.occ=0.6, beta1=0, beta2=0, beta3=0, mean.det=0.3, time.effects=c(0, 0),
  alpha1=0, alpha2=0, alpha3=0, sd.lp=0, b=0))
# psi = 1 (i.e., species occurs at every site)
str(simOcc(mean.occ=1))

# p = 1 (i.e., species is always detected when it occurs)
str(simOcc(mean.det=1))

# Other potentially interesting settings include these:
\donttest{
str(simOcc(J = 2))                 # Only 2 surveys
str(simOcc(M = 1, J = 100))        # No spatial replicates, but 100 measurements
str(simOcc(beta3 = 1))             # Including interaction elev-wind on p
str(simOcc(mean.occ = 0.96))       # A really common species
str(simOcc(mean.occ = 0.05))       # A really rare species
str(simOcc(mean.det = 0.96))       # A really easy species
str(simOcc(mean.det = 0.05))       # A really hard species
str(simOcc(mean.det = 0))          # The dreaded invisible species
str(simOcc(alpha1=-2, beta1=2))    # Opposing effects of elev on psi and p
str(simOcc(J = 10, time.effects = c(-5, 5))) # Huge time effects on p
str(simOcc(sd.lp = 10))            # Huge (random) site effects on p
str(simOcc(J = 10, b = 0))         # No behavioral response in p
str(simOcc(J = 10, b = 2))         # Trap happiness
str(simOcc(J = 10, b = -2))        # Trap shyness
}

# Using categorical covariates only
str(simOccCat())
# Categorical and continuous covariates
str(tmp <- simOccCat(beta1 = -2, beta2 = 2, beta3 = 1,
    mean.detection = 0.3, time.effects = c(-1, 1),
    alpha1 = -1, alpha2 = -3, alpha3 = 0,
    sd.lp = 0.5, b = 2))
# Check how balanced the levels are for HAB
barplot(sort(table(tmp$HAB), decreasing=TRUE), xlab="Habitat category",
    ylab="Frequency", main="Frequency distribution of habitat categories")
}
