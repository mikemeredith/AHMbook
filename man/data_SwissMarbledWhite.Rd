\name{SwissMarbledWhite}
\alias{SwissMarbledWhite}
\encoding{UTF-8}
\docType{data}
\title{
Data for Swiss Marbled White butterflies
}
\description{
Detection/nondetection data for the Marbled White butterfly (\emph{Melanargia galathea}) collected in the canton of Aargau, Switzerland, as part of the Aargau Biodiversity Monitoring Program (LANAG). 519 study plots were each sampled once every five years from 1998 to 2010 with 11 visits each year, resulting in 1337 observed detection histories.
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
Biodiversity Monitoring Program in the canton of Aargau, Switzerland (Courtesy Isabelle Flöss, Abteilung Landschaft und Gewässer des Kanton Aargau). See Appendix C of Roth et al. (2014). %% public domain, CC0

\url{https://figshare.com/articles/dataset/Appendix_C_Field_methods_and_analyses_of_the_butterfly_surveys_of_the_Aargau_Biodiversity_Monitoring_Program_/3559965/1}

}

\references{
Roth, T., Strebel, N. and Amrhein, V. (2014), Estimating unbiased phenological trends by adapting site-occupancy models. \emph{Ecology}, 95: 2144-2154

Kéry, M. & Royle, J.A. (2021) \emph{Applied Hierarchical Modeling in Ecology} AHM2 - 4.13.
}
\seealso{
\code{\link{UKmarbledWhite}}.
}

\examples{
data(SwissMarbledWhite)
str(SwissMarbledWhite)

# Data preparation for AHM2 section 4.13
y <- as.matrix(SwissMarbledWhite[,14:24])    # Grab detection/nondetection data
DATE <- as.matrix(SwissMarbledWhite[,3:13])  # Grab survey dates
for(t in 1:11) {               # Mean-impute date (but don't transform)
  DATE[is.na(DATE[,t]),t] <- mean(DATE[,t], na.rm=TRUE)
}
year <- SwissMarbledWhite$year
nsites <- length(unique(SwissMarbledWhite$site))
nyears <- length(unique(SwissMarbledWhite$year))
nsurveys <- ncol(y)
nobs <- nrow(y)
}
\keyword{datasets}
