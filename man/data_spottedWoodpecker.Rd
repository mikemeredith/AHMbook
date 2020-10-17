\name{spottedWoodpecker}
\alias{spottedWoodpecker}
\encoding{UTF-8}
\docType{data}
\title{
Data for observations of Middle Spotted Woodpecker in Switzerland
}
\description{
The data set is based on checklists which volunteers submit from their bird-watching trips and which were then summarized by day and site (1 km2 quadrat) such that we have the number of surveys (= checklists) per quadrat and day, and the number among these on which a Middle-Spotted Woodpecker was recorded, during 162 days (Julian days 51–212, corresponding to 20 February–31 July). There are data from a total of 144,517 recorded surveys on 116,204 day/quadrat combinations during 26 years from 1545 1 km2 quadrats in which the species was ever recorded in Switzerland since 1985.
}

\usage{data("spottedWoodpecker")}

\format{
\code{spottedWoodpecker} is a data frame with 116,204 rows and 7 columns:

\describe{
  \item{site }{quadrat identifier.}
  \item{coordx, coordy }{the x and y coordinates of the center of the quadrat; the coordinate reference system intentionally not specified.}
  \item{year, jdate }{the year and Julian date of the observations.}
  \item{y }{the number of surveys during which the species was detected.}
  \item{nsurveys }{the number of recorded surveys in the quadrat on that day.}
  }
}

\source{
Database of Swiss Ornithological Institute (courtesy of Nicolas Strebel).
}
\references{
Kéry, M. & Royle, J.A. (2021) \emph{Applied Hierarchical Modeling in Ecology} AHM2 - 4.11.2.
}
\examples{
data(spottedWoodpecker)
str(spottedWoodpecker)
}
\keyword{datasets}
