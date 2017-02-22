\name{SwissSquirrels}
\alias{SwissSquirrels}
\encoding{UTF-8}
\docType{data}
\title{
Data for red squirrels in Switzerland
}
\description{
A file with detection/nondetection data for the red squirrel in 265 1 km2 survey quadrats in Switzerland for 2007, together with covariates. See Examples for code to load the data.
}

% \usage{data("MHB2014")}

\format{
The file \code{SwissSquirrels.txt} is a tab-delimited text file with 265 rows and the following columns:

    \enumerate{
      \item spec.name : the species name.
      \item coordx : the x coordinate of the centre of the quadrat (CRS = CH1903/LV03, EPSG 21781).
      \item coordy : the y coordinate of the centre of the quadrat.
      \item ele : the mean elevation of the quadrat, m.
      \item route : the length of the route walked in the quadrat, km.
      \item forest : percentage forest cover.
      \item det071,	det072,	det073 : 1/0 detection data for 3 survey occasions in 2007.
      \item date071, date072, date073 : Julian date for the 3 survey occasions (1 = 1st April).
      \item dur071, dur072, dur073 : duration of the 3 survey occasions (mins).
    }}

\source{
Data provided with the AHM book at \url{https://www.mbr-pwrc.usgs.gov/pubanalysis/keryroylebook/}
}
\references{
Kéry & Royle (2016) \emph{Applied Hierarachical Modeling in Ecology} Section 10.9
}
\examples{
# To read in the text file:
fn <- file.path(find.package("AHMbook"), "extdata", "SwissSquirrels.txt")
data <- read.table(fn, header = TRUE)
str(data)

}
\keyword{datasets}
