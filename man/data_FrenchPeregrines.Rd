\name{FrenchPeregrines}
\alias{FrenchPeregrines}
\encoding{UTF-8}
\docType{data}
\title{
Data for observations of peregrines from the French Jura
}
\description{
The data are peregrine detection/nondetection data from the French Jura between 1964 and 2016 for 287 cliff sites where a pair had been detected at least once in this period. A large proportion of sites are visited multiple times per year, but unfortunately only the aggregate results are available.
}

\usage{data("FrenchPeregrines")}

\format{
\code{FrenchPeregrines} is a data frame with 287 rows and 56 columns:

\describe{
  \item{site }{cliff identifier.}
  \item{department }{factor, the administrative area (Ain, Jura or Doubs).}
  \item{height }{integer, 1 = low, 2 = medium, 3 = tall cliff.}
  \item{yr1964 to yr2016 }{detection histories for each year: 1 if at least one pair of peregrines was detected, 0 if no pair detected, NA if no survey carried out.}
  }
}

\source{
Groupe Pèlerin Jura (Réne-Jean Monneret)
}
\references{
Kéry & Royle (2019) \emph{Applied Hierarchical Modeling in Ecology} section 15.....
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
