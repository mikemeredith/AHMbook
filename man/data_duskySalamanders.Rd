\name{duskySalamanders}
\alias{duskySalamanders}
\encoding{UTF-8}
\docType{data}
\title{
Age-classified counts of dusky salamanders from two mid-Atlantic states in the USA
}
\description{
Northern dusky salamanders (\emph{Desmognathus fuscus}) were sampled in 12 headwater streams in Maryland and the District of Columbia, USA, annually in June and July from 2005 to 2012. Most streams were sampled at two locations separated by at least 100 m along the stream reach for a total of 21 survey locations. At each location, an observer walked two 15-m transects along the stream-bed and turned cover objects, 6 cm in the water and on the bank within 2 m of the stream. Individuals < 35mm snout-to-vent length were classed as juveniles, including both first-year and second-year juveniles. The data are counts of adults and juveniles.
}

\usage{data("duskySalamanders")}

\format{
\code{duskySalamanders} is four-dimensional (21 x 7 x 2 x 2) array of counts, with dimensions:

\describe{
  \item{site }{a numerical site identifier.}
  \item{year }{the year of the survey, 2005 to 2012.}
  \item{age }{Juveniles or Adults.}
  \item{survey }{June or July.}
}
}
\source{
Evan Grant and colleagues.
}
\references{
Zipkin, E.F. et al (2014) Modeling structured population dynamics using data from unmarked individuals, \emph{Ecology}, 95(1) 22-29.

Kéry, M. & Royle, J.A. (2021) \emph{Applied Hierarchical Modeling in Ecology} AHM2 - 2.9.2
}
\examples{
data(duskySalamanders)
str(duskySalamanders)
}
\keyword{datasets}
