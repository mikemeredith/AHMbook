\name{SwissSquirrels}
\alias{SwissSquirrels}
\encoding{UTF-8}
\docType{data}
\title{
Data for Red Squirrels in Switzerland from the Swiss breeding bird survey MHB
}
\description{
A file with detection/nondetection data for the Red Squirrel (\emph{Sciurus vulgaris}) in 265 1 km2 survey quadrats in Switzerland for 2007, together with covariates (data for the two remaining among the total of 267 quadrats in that monitoring program were not available). See Examples for code to load the data.
}

% \usage{data("MHB2014")}

\format{
The file \code{SwissSquirrels.txt} is a tab-delimited text file with 265 rows and the following columns:

    \itemize{
      \item spec.name : the species name.
      \item coordx : the x coordinate of the center of the quadrat; the coordinate reference system intentionally not specified.
      \item coordy : the y coordinate of the center of the quadrat.
      \item ele : the mean elevation of the quadrat, m.
      \item route : the length of the route walked in the quadrat, km.
      \item forest : percentage forest cover.
      \item det071,	det072,	det073 : 1/0 detection data for 3 survey occasions in 2007.
      \item date071, date072, date073 : Julian date for the 3 survey occasions (1 = 1st April).
      \item dur071, dur072, dur073 : duration of the 3 survey occasions (mins).
    }}

\source{
Swiss Ornithological Institute
}
\references{
Kéry, M. & Royle, J.A. (2016) \emph{Applied Hierarchical Modeling in Ecology} AHM1 - 10.9.
}
\examples{
# To read in the text file
fn <- file.path(find.package("AHMbook"), "extdata", "SwissSquirrels.txt")
data <- read.table(fn, header = TRUE)
str(data)

}
\keyword{datasets}
