\name{waterVoles}
\alias{waterVoles}
\encoding{UTF-8}
\docType{data}
\title{
Data for the Mighty Water Voles of Scotland, UK.
}
\description{
Habitat patches suitable for Water Voles (\emph{Arvicola amphibius}) were surveyed for vole latrines (fecal deposits used for territory marking), used to indicate presence. It is suspected that some of the detections recorded may be false positives. The data set covers 114 patches surveyed up to 3 times per year from 2009 to 2011.
}

\usage{data("waterVoles")}

\format{
\code{waterVoles} is a data frame with 332 rows corresponding to surveys, and 5 columns:

\describe{
  \item{Patch }{an alphanumeric site identifier.}
  \item{y1, y2, y3 }{detection (1)/non-detection (0).}
  \item{Year }{the year of the survey.}
}
}
\source{
Chris Sutherland, Xavier Lambin and colleagues.
}
\references{
Sutherland, C., Elston, D.A., & Lambin, X. (2012) Multi-scale processes in metapopulations: effects of stage structure, rescue effect and correlated extinctions. \emph{Ecology}, 93, 2465-2473.

Sutherland, C.S., Elston, D.A., & Lambin, X. (2014) A demographic, spatially explicit patch occupancy model of metapopulation dynamics and persistence. \emph{Ecology}, 95, 3149-3160.

Kéry, M. & Royle, J.A. (2021) \emph{Applied Hierarchical Modeling in Ecology} AHM2 - 7.2.2.
}
\examples{
data(waterVoles)
str(waterVoles)
}
\keyword{datasets}
