\name{ttdPeregrine}
\alias{ttdPeregrine}
\encoding{UTF-8}
\docType{data}
\title{
Time-to-detection data for peregrines from the French Jura mountains
}
\description{
Between 7 and 9 March 2015, 38 breeding cliffs were visited up to 3 times to survey for Peregrine Falcons (\emph{Falco peregrinus}). Observation duration varied from 3 to 95 minutes. The time to detection was recorded for each bird seen. If no birds were seen, the time was entered as NA.
}

\usage{data("ttdPeregrine")}

\format{
A data frame with 70 rows and the following columns:

\describe{
  \item{Date }{a factor with 3 levels giving the date}
  \item{DayNumber }{the number of the day: 1, 2, or 3}
  \item{SiteNumber }{identification number of the site}
  \item{Start.hour, Start.minute }{the time of starting the search}
  \item{End.hour, End.minute }{the time of ending the search}
  \item{MinOfDay }{the time of the start of the search, minutes after 06:00}
  \item{Tmax }{the duration of the search, in minutes}
  \item{ttd }{the time to detection, minutes; NA if no birds were seen during the search}
  \item{sex }{the sex of the birds seen; NA if no birds were seen}

  }
}

\source{
Based on the field work of one of the authors.
}
\references{
Kéry, M. & Royle, J.A. (2016) \emph{Applied Hierarchical Modeling in Ecology} AHM1 - 10.12.
}
\examples{
# Using the built-in data set instead of the TXT file for the example in Section
#  10.12 of the book is easy, as the format is the same as the file.

?ttdPeregrine  # check the description of the data
data(ttdPeregrine)
# Instead of data <- read.table("ttdPeregrine.txt", header = TRUE, sep = "\t") do:
data <- ttdPeregrine

# Then continue with the rest of the analysis on p.618
}
\keyword{datasets}
