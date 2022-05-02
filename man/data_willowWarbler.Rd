\name{willowWarbler}
\alias{willowWarbler}
\encoding{UTF-8}
\docType{data}
\title{
Constant effort site (CES) data for British Willow Warblers
}
\description{
British Willow Warblers (\emph{Phylloscopus trochilus}) caught at a total of 193 ringing (or banding) stations in England, Wales and Scotland. The willow warbler is a small (8–10 g) migratory species that used to be extremely widespread and common throughout Britain, but whose populations have suffered severe declines since the mid-90s. We analyze data from 11 years (1986–1996) before the major population decline.
}

\usage{data("willowWarbler")}

\format{
\code{willowWarbler} is a list with 4 elements:

\describe{
  \item{\bold{birds }}{a data frame with rows for 10551 birds and the following columns:

    \itemize{
      \item 1986-1996 : capture history, 1 if the bird was captured that year, 0 otherwise.
      \item cesID : numerical code for the Constant Effort Site where the bird was caught.
    }}}
\describe{
  \item{\bold{cells }}{a data frame with rows for 9667 5x5 km cells covering the whole of Britain, and the following columns:

    \itemize{
      \item lon, lat : easting and northing of the center of the cell.
      \item gdd : mean growing degree days (GDD: the sum of daily mean temperatures above 5.5C).
      \item blockID : the ID of the 25x25 km block into which the cell falls.
    }}}
\describe{
  \item{\bold{CES }}{a data frame with rows for each of the Constant Effort Sites, and the following columns:

    \itemize{
      \item cesx, cesy : approximate easting and northing of the site.
      \item blockID : the ID of the 25x25 km block into which the site falls.
      \item cellID : the ID of the 5x5 km cell into which the site falls.
    }}}
\describe{
  \item{\bold{blocks }}{a data frame with rows for 495 25x25 km blocks covering the whole of Britain, and the following columns:

    \itemize{
      \item blockX, blockY : easting and northing of the center of the block.
    }}}
}

\source{
British Ornithological Trust (BTO)
}
\references{
Kéry, M. & Royle, J.A. (2021) \emph{Applied Hierarchical Modeling in Ecology} AHM2 - 3.4.
}
\examples{
data(willowWarbler)
str(willowWarbler)
attach(willowWarbler)
ch <- as.matrix(birds[ , 1:11]) # extract capture histories.

}
\keyword{datasets}
