\name{crossbillAHM}
\alias{crossbillAHM}
\encoding{UTF-8}
\docType{data}
\title{
Data from the Swiss Breeding Bird Survey MHB for European Crossbill from 2001 to 2012
}
\description{
The Swiss breeding bird survey ("Monitoring Häufige Brutvögel" MHB) has monitored the populations of approximately 150 common species since 1999. The MHB sample consists of 267 1-km squares that are laid out as a grid across Switzerland. Fieldwork is conducted by about 200 skilled birdwatchers annually, most of them volunteers. Avian populations are monitored using a simplified territory mapping protocol (also called spot-mapping West of the Atlantic), where each square is surveyed up to three times during the breeding season (only twice above the tree line around 2000 m elevation). Surveys are conducted along a transect that does not change over the years and is 4-6 kms in length on average.

The data frame \code{crossbillAHM} has the data for European Crossbill (\emph{Loxia curvirostra})  from 2001 to 2012.

A variant of this data set for 1999 to 2007 only is included in package \pkg{unmarked}.
}

\usage{data("crossbillAHM")}

\format{
\code{crossbillAHM} is a data frame with 267 rows and 77 columns:

\describe{
  \item{coordx, coordy }{the x and y coordinates of the center of the quadrat; the coordinate reference system intentionally not specified.}
  \item{elev }{the mean elevation of the quadrat, m.}
  \item{forest }{percentage forest cover}
  \item{nsurveys }{the number of replicate surveys \emph{planned} in the quadrat; above the tree-line 2, otherwise 3.}
  \item{det011 to det123 }{1 if the species was detected in the quadrat, 0 otherwise; NA if the corresponding survey was not carried out; the first 2 digits indicate the year and the 3rd digit the survey.}
  \item{date011 to date123 }{Julian date of the survey, 1 = 1st January; the first 2 digits indicate the year and the 3rd digit the survey; NA if the corresponding survey was not carried out.}
  }
}

\source{
Swiss Ornithological Institute
}
\references{
Kéry, M & Royle, J.A. (2021) \emph{Applied Hierarchical Modeling in Ecology} AHM2 - 4.
}
\examples{
data(crossbillAHM)
str(crossbillAHM)

# Extract data as site x survey x year arrays
ytmp <- as.matrix(crossbillAHM[, 6:41]) # matrix of detections
y <- array(ytmp, c(nrow(ytmp), 3, ncol(ytmp)/3))
datetmp <- as.matrix(crossbillAHM[, 42:77])  # matrix of Julian dates
date <- array(datetmp, c(nrow(datetmp), 3, ncol(datetmp)/3))

}
\keyword{datasets}
