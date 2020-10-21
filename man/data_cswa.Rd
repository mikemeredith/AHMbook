\name{cswa}
\alias{cswa}
\encoding{UTF-8}
\docType{data}
\title{
Point count and spot-mapping data for Chestnut-sided Warbler
}
\description{
Chestnut-sided Warblers (\emph{Setophaga pensylvanica}) in 38 \dQuote{wildlife openings} in the White Mountain National Forest were surveyed using point counts and spot-mapping during June and July 2004 as part of Richard Chandler's MS thesis. Counts are for birds detected by song at each point. Two of the larger openings had more than one point, giving 43 points in all. Surveys lasted 10 mins and were divided into intervals of (0–2], (2–5], and (5–10] min. Each point was surveyed 3 times, except for one point surveyed only twice. The spot-mapping data come from intensive surveys of each patch recording individual sightings and nest locations, and using counter-singing to separate adjacent territories.
}

\usage{data("cswa")}

\format{
\code{cwsa} is a list with 3 components:

\describe{
  \item{\bold{counts }}{a points x intervals x occasions array of counts; NAs correspond to the missing survey.}
}

\describe{
  \item{\bold{covs }}{a data frame with rows for 43 points and the following columns:
    \itemize{
      \item time1 : time of start of first survey, hr; NA for the missing survey.
      \item time2 : same for second survey.
      \item time3 : same for third survey.
      \item date1 : date of first survey, NA for the missing survey.
      \item date2 : date of second survey.
      \item date3 : date of third survey.
      \item obs1 : observer ID for first survey.
      \item obs2 : observer ID for second survey.
      \item obs3 : observer ID for third survey.
      \item patchArea : patch area, ha.
      \item plotArea : size of point count area, ha. Most were 50-m radius  but some were smaller because they overlapped adjacent forests which weren’t the focus of the study.
      \item woodHt : woody vegetation height, m.
      \item woodCov : proportion covered by woody vegetation.
    }}}

\describe{
  \item{\bold{spotMaps }}{a data frame with rows for 38 patches and the following columns:
    \itemize{
      \item CSWA : estimated population of chestnut-sided warblers in each patch obtained from spot-mapping.
      \item parea : patch area, ha.
    }}}

}

\source{
David King and colleagues.
}
\references{
Chandler, R.B., King, D.I., & Chandler, C.C. (2009) Effects of management regime on the abundance and nest survival of shrubland birds in wildlife openings in northern New England, USA. \emph{Forest Ecology and Management}, 258, 1669-1676.

Chandler, R.B., Royle, J.A., & King, D.I. (2011) Inference about density and temporary emigration in unmarked populations. \emph{Ecology}, 92, 1429-1435.

Kéry, M. & Royle, J.A. (2021) \emph{Applied Hierarchical Modeling in Ecology} AHM2 - 2.
}
\examples{
data(cswa)
str(cswa)

# Replace the NAs in the covariates with the mean from other surveys:
covs <- cswa$covs
covs$time1[is.na(covs$time1)] <- mean(covs$time1, na.rm=TRUE)
covs$date1[is.na(covs$date1)] <- mean(covs$date1, na.rm=TRUE)
# Now use covs instead of cswa$covs

}
\keyword{datasets}
