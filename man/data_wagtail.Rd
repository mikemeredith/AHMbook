\name{wagtail}
\alias{wagtail}
\encoding{UTF-8}
\docType{data}
\title{
Data for Dutch wagtails
}
\description{
The Dutch Centre for Field Ornithology (Sovon) %% "Centre" is correct
monitored grassland birds in Flevoland between April and mid-July 2011. 235 points were surveyed on up to 4 occasions; observations were divided into distance classes 50m wide, and the number of observations in each class recorded. These data are for the Yellow Wagtail (\emph{Motacilla flava}).
}

\usage{data("wagtail")}

\format{
\code{wagtail} is a list with 8 elements:

\describe{
  \item{potato }{for each point, the percentage of the area which are potato fields.}
  \item{grass }{for each point, the percentage of the area under permanent grassland.}
  \item{lscale }{for each point, an index of whether the landscape is open (0) or closed (100).}
  \item{date }{a points x occasions matrix, the Julian date of the survey.}
  \item{hour }{a points x occasions matrix, the hour of the survey.}
  \item{breaks }{the boundaries between the distance classes; birds more than 300m from the point were not included.}
  \item{Y }{a matrix of counts, with a row for each site; columns 1 to 6 give the counts in the distance classes for the 1st survey occasion, columns 7 to 12 for the 2nd occasion, and so on (this is the format required for \code{unmarkedFrameGDS}).}
  \item{rep }{a points x occasions character matrix with the occasion number (this is used as a categorical variable in the analysis).}
}
}
\source{
Dutch Centre for Field Ornithology (SOVON) %% 'Centre' is correct
}
\references{
Kéry, M. & Royle, J.A. (2016) \emph{Applied Hierarchical Modeling in Ecology} AHM1 - 9.5.3.
}
\examples{
data(wagtail)
str(wagtail)
}
\keyword{datasets}
