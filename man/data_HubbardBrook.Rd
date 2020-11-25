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
  \item{\bold{counts }}{a sites x replicates x years x species array of counts, the number of unique individuals detected within 50m of the point count location.}

  \item{\bold{sitecov }}{a data frame with rows for 373 sites and the following columns:

    \itemize{
      \item PlotID : a numeric site identifier.
      \item Longitude : the longitude of the point (WGS84); 2 sites have missing data.
      \item Latitude : the latitude of the point (WGS84); 2 sites have missing data.
      \item Elev : the elevation of the point, m.
      \item Aspect : aspect of the sample location (degrees).
      \item Slope : slope of the sample location.
      \item Forest : general forest cover type of sample location.
    }}

  \item{\bold{dates }}{a sites x replicates x years array of ordinal day of count (1 Jan = 1); NA entries occur when surveys were not carried out.}

  \item{\bold{times }}{a sites x replicates x years array with the start time of the survey, hours after midnight.}
}
}

\source{
Rodenhouse N.L. & Sillett, T.S. (2019) Valleywide Bird Survey, Hubbard Brook Experimental Forest, 1999-2016 (ongoing). Environmental Data Initiative. <https://doi.org/10.6073/pasta/faca2b2cf2db9d415c39b695cc7fc217>. Dataset accessed 2020-01-07.}

\references{
Betts, M.G., Rodenhouse, N.L., Sillett, T.S., Doran, P.J. & Holmes, R.T. (2008). Dynamic occupancy models reveal within-breeding season movement up a habitat quality gradient by a migratory songbird. \emph{Ecography} 31:592–600.

Goetz, S.J., Steinberg, D., Betts, M.G., Holmes, R.T., Doran, P.J., Dubayah, R., & Hofton, M. (2010). Lidar remote sensing variables predict breeding habitat of a Neotropical migrant bird. \emph{Ecology} 91:1569–1576.

Van Tatenhove, A., Filiberti, E., Sillett, T.S., Rodenhouse, N.L. & Hallworth, M. T. (2019). Climate-related distribution shifts of migratory songbirds and sciurids in the White Mountain National Forest. \emph{Forests} 10:84.

Kéry, M. & Royle, J.A. (2021) \emph{Applied Hierarchical Modeling in Ecology} AHM2 - 8.2.
}

\examples{
data(HubbardBrook)
str(HubbardBrook)
}
\keyword{datasets}
