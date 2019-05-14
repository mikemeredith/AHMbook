\name{HubbardBrook}
\alias{HubbardBrook}
\encoding{UTF-8}
\docType{data}
\title{
Warbler point count data from the Hubbard Brook Experimental Forest
}
\description{
Three replicated point counts of 10 min each, out to a maximum distance of 100 m, were conducted in each spring between 1999 and 2018 at a total of 373 locations in the Hubbard Brook Experimental Forest in New Hampshire. This data set contains the data for the following 13 species:
\itemize{
  \item AMRE - American Redstart
  \item BAWW - Black-and-white Warbler
  \item BHVI - Blue headed Vireo
  \item BLBW - Blackburnian Warbler
  \item BLPW - Blackpoll Warbler
  \item BTBW - Black-throated Blue Warbler
  \item BTNW - Black-throated Green Warbler
  \item CAWA - Canada Warbler
  \item MAWA - Magnolia Warbler
  \item MYWA - Myrtle Warbler
  \item NAWA - Nashville Warbler
  \item OVEN - Ovenbird
  \item REVI - Red-eyed Vireo
}
}

\usage{data("HubbardBrook")}

\format{
\code{HubbardBrook} is a list with 4 elements:

\describe{
  \item{\bold{counts }}{a sites x replicates x years x species array of counts, the number of unique individuals counted within 50m of the point count location.}

  \item{\bold{sitecov }}{a data frame with rows for 373 sites and the following columns:

    \itemize{
      \item PlotID : a numeric site identifier.
      \item Longitude : the longitude of the point (WGS84); 2 sites have missing data.
      \item Latitude : the latitude of the point (WGS84); 2 sites have missing data.
      \item Elev : the elevation of the point, m.
      \item Aspect : aspect of the sample location (radians).
      \item Slope : slope of the sample location.
      \item Forest : general forest cover type of sample location.
    }}

  \item{\bold{dates }}{a sites x replicates x years array of ordinal day of count (1 Jan = 1); NA entries occur when surveys were not carried out.}

  \item{\bold{times }}{a sites x replicates x years array with the start time of the survey, hours after midnight.}
}
}

\source{
Scott Sillett ... insert proper attribution.
}
\references{
Zipkin...??

Kéry & Royle (2019) \emph{Applied Hierarchical Modeling in Ecology} Section 19.2
}
\examples{
data(HubbardBrook)
str(HubbardBrook)
}
\keyword{datasets}
