\name{BerneseOberland}
\alias{BerneseOberland}
\docType{data}
\title{
Landscape data for the Bernese Oberland around Interlaken, in the Alps of central Switzerland
}
\description{
Spatially-referenced data on elevation, forest cover, and water cover at a 1km-sq resolution for a 50km x 50km region in Switzerland, the Bernese Oberland. This is a small subset of the \code{Switzerland} data in the \pkg{unmarked} package. See the Examples for a plot of the location.
}
\usage{data(BerneseOberland)}
\format{
  A data frame with 2500 observations on the following 5 variables.
  \describe{
    \item{x }{Easting (m)}
    \item{y }{Northing (m)}
    \item{elevation }{a numeric vector (m)}
    \item{forest }{a numeric vector (percent cover)}
    \item{water }{a numeric vector (percent cover)}
  }
}
\details{
Forest and water coverage (in percent area) was computed using the
1992-97 landcover dataset of the Swiss Federal Statistical Office
(http://www.bfs.admin.ch). Median elevation (in meters) was computed
using a median aggregation of the digital elevation model of the Swiss
Federal Statistical Office.

x and y are the coordinates of the center of each 1km2 pixel.

The coordinate reference system is intentionally not specified.

These data can only be used for non-profit projects. Otherwise, written
permission must be obtained from the Swiss Federal Statistical Office.
}
\source{
Swiss Federal Statistical Office (http://www.bfs.admin.ch)
}
\examples{
data(BerneseOberland)
str(BerneseOberland)

\donttest{
library(raster)
data(Switzerland, package="unmarked")
r1 <- with(Switzerland, rasterFromXYZ(data.frame(x = x, y = y, z = elevation)))
mapPalette1 <- colorRampPalette(c("grey", "yellow", "orange", "red"))
plot(r1, col = mapPalette1(100), axes = FALSE, box = FALSE, zlim = c(0, 4500),
  main = "Location of Bernese Oberland data set")
with(BerneseOberland, rect(min(x), min(y), max(x), max(y), lwd=2))
}

}
\keyword{datasets}
