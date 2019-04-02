\name{SwissEagleOwls}
\alias{SwissEagleOwls}
\encoding{UTF-8}
\docType{data}
\title{
Data for Eagle Owls in Switzerland
}
\description{
Data stem from more or less opportunistic surveys of a total of 274 territories of the Eagle Owl (\emph{Bubo bubo}), the largest owl in the world, between 2007 and 2016 (10 years). Sites are all over Switzerland and the complete list of sites where a territorial Eagle Owl was ever detected during this period. For site/year combinations with more than 20 observations, 20 were randomly selected for inclusion.
}

\usage{data("SwissEagleOwls")}

\format{
\code{SwissEagleOwls} is a list with 2 elements:

\describe{
  \item{\bold{obs }}{a data frame with a row for each of the 5974 observations and the following columns:

    \enumerate{
      \item site_name : a numeric site ID.
      \item year : the year of the observation.
      \item jdate : the Julian date of the observation (1 = Jan 1st).
      \item y : the observed state: 0 = species not detected; 1 = detection of a single bird; 2 = detection of a pair without evidence of breeding; 3 = detection of a pair with evidence of breeding.
    }}}
\describe{
  \item{\bold{sites }}{a data frame with rows for 274 1x1km quadrats containing the sites, and the following columns:

    \enumerate{
      \item site_name : the numeric site ID, corresponding to site_name in the obs data frame.
      \item region : codes for six regions of Switzerland.
      \item elev : the elevation of the centre of the quadrat, m.
      \item forest : proportion of the quadrat covered by forest.
      \item slope : the slope of the quadrat in degrees.
      \item roads : total road length in the quadrat, m.
    }}}
}

\source{
Swiss Ornithological Institute
}
\references{
Kéry & Royle (2016) \emph{Applied Hierarchical Modeling in Ecology} Chapter 17.
}
\examples{
data("SwissEagleOwls")
str(SwissEagleOwls)
}
\keyword{datasets}
