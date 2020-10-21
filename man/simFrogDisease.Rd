\name{simFrogDisease}
\alias{simFrogDisease}
\encoding{UTF-8}
\title{
Simulate detection data for diseased frogs
}
\description{
Simulates detection data for an amphibian species affected by the fungal pathogen \emph{Batrachochytrium dendrobatidis} (Lips et al., 2006). Surveys are carried out on  multiple occasions each year. Observers count frogs and record whether they were infected or not. Over time, individuals are subject to mortality and there may also be recruitment of new individuals into the population. Individuals transition between the two states: infected frogs may become uninfected (i.e., can shed the fungus) and vice versa. It also seems reasonable to expect that survival probability differs between the two states.
}
\usage{
simFrogDisease(nsites = 100, nyears = 3, nsurveys = 3, alpha.lam = 3,
    omega = c(0.9, 0.7), gamma = c(2,1), p = c(0.8, 0.8, 0.8),
    recovery = 0.1, infection = 0.1)

}
\arguments{
  \item{nsites}{
the number of sites surveyed.
}
  \item{nyears}{
the number of years of the study.
}
  \item{nsurveys}{
the number of surveys each year.
}
  \item{alpha.lam}{
mean abundance per site in the first year.
}
  \item{omega}{
vector length 2, state-specific survival, noninfected and infected.
}
  \item{gamma}{
vector length 2, state-specific recruitment, noninfected and infected.
}
  \item{p}{
vector length \code{nyears}, probability of detection of individuals.
}
  \item{recovery}{
probability of recovery for an infected frog.
}
  \item{infection}{
probability of infection for a noninfected frog.
}
}
\value{
A list with the values of the arguments and the following additional elements:
\item{SN }{sites x intervals, number of noninfected frogs surviving}
\item{SI }{sites x intervals, number of infected frogs surviving}
\item{GN }{sites x intervals, number of noninfected frogs recruited}
\item{GI }{sites x intervals, number of infected frogs recruited}
\item{TrI }{sites x intervals, number of infected frogs recovering}
\item{TrN }{sites x intervals, number of noninfected frogs becoming infected}
\item{NN }{sites x years, number of noninfected frogs at each site}
\item{NI }{sites x years, number of infected frogs at each site}
\item{yN }{sites x years x surveys, number of noninfected frogs detected}
\item{yI }{sites x years x surveys, number of infected frogs detected}
}
\references{
Kéry, M. & Royle, J.A. (2021) \emph{Applied Hierarchical Modeling in Ecology} AHM2 - 2.9.1.

Lips, K.R., Brem, F., Brenes, R., Reeve, J.D., Alford, R.A., Voyles, J., et al., 2006. Emerging infectious disease and the loss of biodiversity in a Neotropical amphibian community. \emph{Proc. Natl. Acad. Sci. USA}, 103, 3165-3170.
}
\author{
Marc Kéry and Andy Royle
}

\examples{
# Generate a simulated data set with default arguments and look at the structure:
tmp <- simFrogDisease()
str(tmp)
}
