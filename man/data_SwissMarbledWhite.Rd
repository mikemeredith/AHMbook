\name{SwissMarbledWhite}
\alias{SwissMarbledWhite}
\encoding{UTF-8}
\docType{data}
\title{
Data for Swiss Marbled White butterflies
}
\description{
Detection/nondetection data for the Marbled White butterfly (\emph{Melanargia galathea}) collected from 1998 to 2010 at a total of 519 study plots in Switzerland, as part of the Swiss Biodiversity Monitoring, which were each sampled once every five years with 11 visits each year, resulting in 1337 observed detection histories.
}

\usage{data("SwissMarbledWhite")}

\format{
A data frame with 1337 rows and the following columns:

\describe{
  \item{site }{identification number for site}
  \item{year }{the year of the record}
   \item{Day1, Day2, Day3, Day4, Day5, Day6, Day7, Day8, Day9, Day10, Day11 }{Julian date of the survey, 1 = 1st April. NA indicates that the date of the survey was not recorded; see Examples for a way to impute dates.}
   \item{y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11 }{1 if the species was recorded during the survey, 0 otherwise.}
  }
}

\source{
Aargau Biodiversity Monitoring Program (LANAG) (Courtesy Tobias Roth). See Appendix C of Roth et al (2014).
}
\references{
Roth, T., Strebel, N. and Amrhein, V. (2014), Estimating unbiased phenological trends by adapting site-occupancy models. \emph{Ecology}, 95: 2144-2154

Kéry & Royle (2019) \emph{Applied Hierarchical Modeling in Ecology} Section 15...
}
\seealso{
\code{\link{UKmarbledWhite}}.
}

\examples{
data(SwissMarbledWhite)
str(SwissMarbledWhite)

# Impute values for the missing dates:
smw <- SwissMarbledWhite
for(t in 3:13)
  smw[is.na(smw[, t]), t] <- mean(smw[, t], na.rm=TRUE)
summary(smw)
}
\keyword{datasets}
