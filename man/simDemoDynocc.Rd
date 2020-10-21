\name{simDemoDynocc}
\alias{simDemoDynocc}
\encoding{UTF-8}
\title{
Simulate data under a demographic dynamic occupancy model
}
\description{
Function to simulate detection/nondetection data under a variant of the demographic occupancy (or 'local survival') model of Roth & Amrhein (2010). Data are simulated in an 'unconditional' manner, i.e., for each site from first to last year. All parameters can be made year-dependent by specification of a range within which annual values will be drawn from uniform distributions.
}
\usage{
simDemoDynocc(nsites = 100, nyears = 10, nvisits = 5, psi1 = 0.6,
    range.phi = c(0.2, 0.9), range.r = c(0, 0.4), range.p = c(0.1, 0.9),
    show.plot=TRUE)
}
\arguments{
  \item{nsites}{
Number of sites.
}
  \item{nyears}{
Number of years (or 'seasons', as they are somewhat confusingly often called in the occupancy literature).
}
  \item{nvisits}{
Number of replicate surveys (= occasions) within a year.
}
  \item{psi1}{
occupancy probability in first year.
}
  \item{range.phi}{
bounds of uniform distribution from which annual local probability of persistence is randomly drawn.
}
  \item{range.r}{
bounds of uniform distribution from which annual local probability of colonization is randomly drawn.
}
  \item{range.p}{
bounds of uniform distribution from which annual probability of detection is randomly drawn.
}
  \item{show.plot}{
If TRUE, plots of results are displayed; set to FALSE if running simulations.
}
}

\value{
A list with the values of the arguments input and the following additional elements:
  \item{phi }{persistence for each interval, a vector of length nyears - 1}
  \item{r }{colonization for each interval, a vector of length nyears - 1}
  \item{p }{detection probability for each year, a vector of length nyears}
  \item{z }{true occurrence state, a nsites x nyears matrix of 0/1}
  \item{y }{the observed detection history, a nsites x nvisits x nyears array}
  \item{f }{year of first detection, a vector of length nsites}
  \item{nocc.true }{the true number of occupied sites, a vector of length nyears}
  \item{nocc.true }{the observed number of occupied sites, a vector of length nyears}
}
\references{
Roth, T. & Amrhein, V. (2010), Estimating individual survival using territory occupancy data on unmarked animals. \emph{Journal of Applied Ecology}, 47, 386-392.

Kéry, M. & Royle, J.A. (2021) \emph{Applied Hierarchical Modeling in Ecology} AHM2 - 4.
}
\author{
Marc Kéry & Andy Royle
}

\examples{
# Generate data with the default arguments and look at the structure:
str(data <- simDemoDynocc() )                # Implicit defaults
str(data <- simDemoDynocc(psi1 = 1))         # All sites initially occupied
str(data <- simDemoDynocc(nsites = 1000))    # Plenty more sites
str(data <- simDemoDynocc(nyears = 100))     # Plenty more years
str(data <- simDemoDynocc(nvisits = 20))     # Plenty more visits
str(data <- simDemoDynocc(range.phi = c(0.8, 0.8))) # Constant survival
str(data <- simDemoDynocc(range.phi = c(0.2,0.3), range.r = c(0,0.2))) # Decline
str(data <- simDemoDynocc(range.phi = c(0.8,1), range.r = c(0.5,0.7))) # Increase
str(data <- simDemoDynocc(nvisits = 1))      # Single visit
str(data <- simDemoDynocc(range.p = c(1,1))) # Perfect detection
}

