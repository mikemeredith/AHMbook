\name{SwissEagleOwls}
\alias{SwissEagleOwls}
\encoding{UTF-8}
\docType{data}
\title{
Territory-level, multi-state detection/nondetection data for Eagle Owls in Switzerland
}
\description{
Multi-state detection/nondetection data stem from more or less opportunistic surveys in a total of 274 territories of the Eurasian Eagle Owl (\emph{Bubo bubo}), the largest owl in the world, between 2007 and 2016 (10 years). Sites are all over Switzerland and represent the complete list of sites where a territorial Eagle Owl was ever detected during this period. For site/year combinations with more than 20 observations, 20 were randomly selected for inclusion.

The data set recognizes four observed states: 'nondetection', 'detection of single bird', 'detection of a pair without breeding evidence' and 'detection of a pair with confirmed breeding'; these being the classical criteria of bird atlases for possible, probable or confirmed breeding. For the analysis in the book, the last two are aggregated into a single state 'Pair'.
}

\usage{data("SwissEagleOwls")}

\format{
\code{SwissEagleOwls} is a list with 2 elements:

\describe{
  \item{\bold{obs }}{a data frame with a row for each of the 5974 observations and the following columns:

    \itemize{
      \item site_name : a numeric site ID.
      \item year : the year of the observation.
      \item jdate : the Julian date of the observation (1 = Jan 1st).
      \item y : the observed state: 0 = species not detected; 1 = detection of a single bird; 2 = detection of a pair without evidence of breeding; 3 = detection of a pair with evidence of breeding.
    }}}
\describe{
  \item{\bold{sites }}{a data frame with rows for 274 1x1km quadrats containing the sites, and the following columns:

    \itemize{
      \item site_name : the numeric site ID, corresponding to site_name in the obs data frame.
      \item region : codes for six regions of Switzerland.
      \item elev : the elevation of the center of the quadrat, m.
      \item forest : proportion of the quadrat covered by forest.
      \item slope : the slope of the quadrat in degrees.
      \item roads : total road length in the quadrat, m.
    }}}
}

\source{
Swiss Ornithological Institute
}
\references{
Kéry, M. & Royle, J.A. (2021) \emph{Applied Hierarchical Modeling in Ecology} AHM2 - 6.
}
\examples{
data("SwissEagleOwls")
str(SwissEagleOwls)
}
\keyword{datasets}
