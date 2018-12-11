\name{middleSpottedWoodpecker}
\alias{middleSpottedWoodpecker}
\encoding{UTF-8}
\docType{data}
\title{
Data for observations of middle spotted woodpecker in Switzerland
}
\description{
Marc is preparing a description to put here.
}

\usage{data("middleSpottedWoodpecker")}

\format{
\code{middleSpottedWoodpecker} is a data frame with 116204 rows and 7 columns:

\describe{
  \item{site }{quadrat identifier.}
  \item{coordx, coordy }{the x and y coordinates of the centre of the quadrat; the coordinate reference system intentionally not specified.}
  \item{year, jdate }{the year and Julian date of the observations.}
  \item{y }{the number of surveys during which the species was detected.}
  \item{nsurveys }{the number of recorded surveys in the quadrat on that day.}
  }
}
\details{

May not need this section.
}

\source{
Database of Swiss Ornithological Institute (courtesy of Nicolas Strebel).
}
\references{
Kéry & Royle (2019) \emph{Applied Hierarchical Modeling in Ecology} chapter 15.
}
\examples{
data(middleSpottedWoodpecker)
str(middleSpottedWoodpecker)
}
\keyword{datasets}
