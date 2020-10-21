\name{treeSparrow}
\alias{treeSparrow}
\encoding{UTF-8}
\docType{data}
\title{
American Tree Sparrow data from Alaska
}
\description{
Mizel et al. (2018) describe the data on American Tree Sparrows (\emph{Spizella arborea}) that were collected in May and June 2016 at Noatak National Preserve, Alaska. The study area was a 10.8 x 9.6 km rectangle divided up into 288 plots of 600 m x 600 m, of which 150 had some coverage by a transect running through them. Routes were usually chosen to intersect the centroid of the cell and deflected by 0 or 90 degrees when required to move into an adjacent cell. In some cases, a transect was shifted to avoid terrain features. Distance measurements were truncated at 250m from the observer. Each visited plot was surveyed an average of three times (max = 5). Each plot was gridded to a 30-m resolution and spatial attributes were associated with each of these pixels. This resulted in an average of 328 pixels per plot within the buffered transect. Bird locations were assigned to the nearest pixel centroid.
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
      \item juldate : the Julian date of the survey, day 1 = 1 Jan.
      \item effort : time spent in each plot during each visit (in hours)
      \item reltime : number of hours since 0200 (there is no sunrise during the survey period in the Arctic)
      \item Noise : background noise (combination of wind, creek, and mosquito-associated noise): 0 = no noise, 1 = slight, 2 = noticeable reduction in hearing, 3 = prohibitive noise
      \item Observer2 : Observed ID
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
      \item NEAR_DIST : distance between the center of the pixel and the transect line
      \item NDVI : Normalized Difference Vegetation Index
      \item elev : elevation of the center of the pixel, m.
    }}
}
}

\source{
Mizel et al. (2018)
}
\references{
Mizel, J.D., Schmidt, J.H. and Lindberg, M.S. (2018) Accommodating temporary emigration in spatial distance sampling models. \emph{Journal of Applied Ecology}, 55, 1456-1464.

Kéry, M. & Royle, J.A. (2021) \emph{Applied Hierarchical Modeling in Ecology} AHM2 - 11.8.4.
}
\examples{
data(treeSparrow)
str(treeSparrow)
}
\keyword{datasets}
