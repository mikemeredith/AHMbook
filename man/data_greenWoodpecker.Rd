\name{greenWoodpecker}
\alias{greenWoodpecker}
\encoding{UTF-8}
\docType{data}
\title{
Count data from the Swiss Breeding Bird Survey for Green Woodpeckers from 2004 to 2017
}
\description{
The Swiss breeding bird survey ("Monitoring Häufige Brutvögel" MHB) has monitored the populations of approximately 150 common species since 1999. The MHB sample consists of 267 1-km squares that are laid out as a grid across Switzerland. Fieldwork is conducted by about 200 skilled birdwatchers annually, most of them volunteers. Avian populations are monitored using a simplified territory mapping protocol (also called spot-mapping West of the Atlantic), where each square is surveyed up to three times during the breeding season (only twice above the tree line around 2000 m elevation). Surveys are conducted along a transect that does not change over the years and is 4-6 kms in length on average.

The data frame \code{greenWoodpecker} has count data for Green Woodpeckers \emph{Picus viridis} from 2004 to 2017. There are some missing values: see Details.
}

\usage{data("greenWoodpecker")}

\format{
\code{greenWoodpecker} is a data frame with 267 rows and 146 columns:

\describe{
  \item{x, y }{the x and y coordinates of the center of the quadrat; the coordinate reference system intentionally not specified.}
  \item{nsurveys }{the number of replicate surveys \emph{planned} in the quadrat; above the tree-line 2, otherwise 3.}
  \item{elevation }{the mean elevation of the quadrat, m.}
  \item{forest }{percentage forest cover}
  \item{route.length }{the length of the route}
  \item{X2004.1.count to X2017.3.count }{number of territories observed for each survey; 2004.1 indicates the first survey in 2004.}
  \item{X2004.1.date to X2017.3.date }{Julian date of the survey, 1 = a weekend near 15 April.}
  \item{X2004.1.time to X2017.3.time }{time spent in the quadrat during the survey, mins.}
  \item{X2004.observer to X2017.observer }{numerical code identifying the observer.}
  }
}
\details{

Missing values in the \code{count} and \code{date} columns indicate that the corresponding survey was not carried out.

On 26 occasions when surveys were carried out, the time taken was not recorded, resulting in additional NAs in the \code{time} columns.
}

\source{
Swiss Ornithological Institute
}
\references{
Kéry, M. & Royle, J.A.  (2021) \emph{Applied Hierarchical Modeling in Ecology} AHM2 - 2.2.
}
\examples{
data(greenWoodpecker)
str(greenWoodpecker)

}
\keyword{datasets}
