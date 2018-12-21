\name{simDemoDynocc}
\alias{simDemoDynocc}
\encoding{UTF-8}
\title{
Simulate data under a demographic dynamic occupancy model
}
\description{
Function to simulate detection/nondetection data under a variant of the demographic occupancy (or 'local survival') model of Roth & Amrhein (2010). Data are simulated in an 'unconditional' manner, i.e., stochastic  modeling starts after a territory is first detected. All parameters can be made year-dependent by specification of a range within which annual values will be drawn from uniform distributions.
}
\usage{
simDemoDynocc(nsite = 100, nyear = 10, nvisit = 5, psi1 = 0.6,
    range.phi = c(0.2, 0.9), range.r = c(0, 0.4), range.p = c(0.1, 0.9),
    show.plot=TRUE)
}
\arguments{
  \item{nsite}{
Number of sites.
}
  \item{nyear}{
Number of years (or 'seasons').
}
  \item{nvisit}{
Number of replicate surveys (= occasions) within a year.
}
  \item{psi1}{
occupancy probability in first year.
}
  \item{range.phi}{
bounds of uniform distribution from which annual local probability of persistence is randomly drawn.
}
  \item{range.r}{
bounds of uniform distribution from which annual local probability of colonisation is randomly drawn.
}
  \item{range.p}{
bounds of uniform distribution from which annual probability of detection is randomly drawn.
}
  \item{show.plot}{
If TRUE, plots of results are displayed; set to FALSE if running many simulations.
}
}

\value{
A list with the values of the arguments input and the following additional elements:
  \item{phi }{persistence for each interval, a vector of length nyear - 1}
  \item{r }{colonisation for each interval, a vector of length nyear - 1}
  \item{p }{detection probability for each year, a vector of length nyear}
  \item{z }{true occurrence state, a nsite x nyear matrix of 0/1}
  \item{y }{the observed detection history, a nsite x nvisit x nyear array}
  \item{f }{year of first detection, a vector of length nsite}
  \item{nocc.true }{the true number of occupied sites, a vector of length nyear}
  \item{nocc.true }{the observed number of occupied sites, a vector of length nyear}
}
\references{
Roth, T. and Amrhein, V. (2010), Estimating individual survival using territory occupancy data on unmarked animals. \emph{Journal of Applied Ecology}, 47: 386-392.

Kéry & Royle (2017) \emph{Applied Hierarchical Modeling in Ecology} Chapter 15
}
\author{
Marc Kéry
}

\examples{
# Generate data with the default arguments and look at the structure:
str(data <- simDemoDynocc() )            # Implicit defaults
str(data <- simDemoDynocc(psi1 = 1))     # All sites initially occupied
str(data <- simDemoDynocc(nsite = 1000)) # Plenty more sites
str(data <- simDemoDynocc(nyear = 100))  # Plenty more years
str(data <- simDemoDynocc(nvisit = 20))  # Plenty more visits
str(data <- simDemoDynocc(range.phi = c(0.8, 0.8))) # Constant survival
str(data <- simDemoDynocc(range.phi = c(0.2,0.3), range.r = c(0,0.2))) # Decline
str(data <- simDemoDynocc(range.phi = c(0.8,1), range.r = c(0.5,0.7))) # Increase
str(data <- simDemoDynocc(nvisit = 1))   # Single visit
str(data <- simDemoDynocc(range.p = c(1,1))) # Perfect detection
}

