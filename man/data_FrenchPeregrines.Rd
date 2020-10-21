\name{FrenchPeregrines}
\alias{FrenchPeregrines}
\encoding{UTF-8}
\docType{data}
\title{
Data for observations of peregrines from the French Jura mountains
}
\description{
The data are detection/nondetection data of the Peregrine Falcon (\emph{Falco peregrinus}) from the wild and wonderful French Jura between 1964 and 2016 for 284 cliff sites (= territories, or sites in the context of a site-occupancy model) where a pair had been detected at least once in this period. A large proportion of sites are visited multiple times per year, but unfortunately only the aggregate results are available in each year, i.e., whether a pair was detected at least once, or never.
}

\usage{data("FrenchPeregrines")}

\format{
\code{FrenchPeregrines} is a data frame with 284 rows and 56 columns:

\describe{
  \item{site }{cliff (or site) identifier.}
  \item{department }{factor, the administrative area (Ain, Jura or Doubs).}
  \item{height }{factor, height of the cliff, low, medium, or tall.}
  \item{yr1964 to yr2016 }{detection histories for each year: 1 if a pair of peregrines was detected during at least one survey, 0 if no pair was detected, NA if no survey was carried out in that year.}
  }
}

\source{
Groupe Pèlerin Jura (René-Jean Monneret, René Ruffinoni, and colleagues)
}
\references{
Kéry, M. & Royle, J.A. (2021) \emph{Applied Hierarchical Modeling in Ecology} AHM2 - 4.11.
}
\examples{
data(FrenchPeregrines)
str(FrenchPeregrines)
# Extract the capture history data:
ch <- as.matrix(FrenchPeregrines[, 4:56])
dim(ch)
range(ch, na.rm=TRUE)
}
\keyword{datasets}
