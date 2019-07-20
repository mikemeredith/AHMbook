\name{treeSparrow}
\alias{treeSparrow}
\encoding{UTF-8}
\docType{data}
\title{
Tree sparrow data
}
\description{
describe data collection
}

\usage{data("treeSparrow")}

\format{
\code{treeSparrow} is a list with 3 elements:

\describe{

  \item{\bold{surveyData }}{a data frame with rows for 466 surveys and the following columns:

    \itemize{
      \item Site : a numeric site identifier, 1 to 150
      \item Visit : a numeric visit identifier, 1 to 5
      \item Survey : a numeric survey identifier, 1 to 466
      \item count : the number of tree sparrows detected.
      \item juldate : the Julian date of the survey. ?? origin??
      \item effort : ??
      \item reltime : ??
      \item Noise : ??
      \item Observer2 : ??
    }}

  \item{\bold{obs }}{a data frame with rows for 325 birds observed and the following columns:

    \itemize{
      \item Pixel : a numeric pixel identifier, 1 to 49250
      \item Site : a numeric site identifier, 1 to 150
      \item Visit : a numeric visit identifier, 1 to 5
      \item NEAR_DIST :  distance between the bird and the transect line
      \item SurveyID : a numeric survey identifier, 1 to 466
      \item ind : an index identifying the individual within the survey.
    }}

  \item{\bold{pixels }}{a data frame with rows for 49250 pixels and the following columns:

    \itemize{
      \item pixID : a numeric pixel identifier, 1 to 49250
      \item Site : a numeric site identifier, 1 to 150
      \item X : x-coordinate of the pixel
      \item Y : y-coordinate of the pixel
      \item NEAR_DIST : distance between the centre of the pixel and the transect line
      \item NDVI : Normalized Difference Vegetation Index
      \item elev : elevation of the centre of the pixel, m.
    }}
}
}

\source{
Mizel et al. (2018) ...?? insert proper attribution.
}
\references{
Mizel et al. (2018)...??

Kéry & Royle (2019) \emph{Applied Hierarchical Modeling in Ecology} Section 11.8.4
}
\examples{
data(treeSparrow)
str(treeSparrow)
}
\keyword{datasets}
