\name{Color Ramps}
\alias{Color_Ramps}
\alias{rampYOR}
\alias{rampGreens}
\alias{rampGreys}
\alias{rampBYR}
\alias{rampGBr}
\encoding{UTF-8}
\title{
Color Ramps
}
\description{
These create vectors of \code{n} contiguous colors for use in image and raster plots. They are similar to heat.colors or terrain.colors but with extra arguments to give greater control of colors. Currently included ramps are:

* sequential yellow-orange-red

* sequential light to dark green

* sequential light to dark gray

* diverging blue-yellow-red

* diverging green-brown

Sequential ramps go from light to dark, diverging ramps from dark through light to dark. All these ramps are color-blind safe and based on color palettes from \url{https://colorbrewer2.org/} by Cynthia Brewer and colleagues.
}
\usage{
rampYOR(n=5, range=1:9, bias=1, ...)
rampGreens(n=5, range=1:9, bias=1, ...)
rampGreys(n=5, range=1:9, bias=1, ...)
rampBYR(n=5, range=1:11, bias=1, ...)
rampGBr(n=5, range=1:11, bias=1, ...)

}
\arguments{
  \item{n}{
the number of colors to be in the palette.
}
  \item{range}{
the range of colors from the underlying palette to include: each ramp is based on a palette of 9 (if sequential) or 11 (if diverging) colors; range allows the extreme colors to be excluded.
}
  \item{bias}{
a positive number: values >1 give more widely spaced colors at the high end.
}
  \item{\dots}{
additional arguments passed to \code{\link{colorRampPalette}}.
}
}

\value{
Return a vector of color values in hexadecimal format.
}
\author{
Mike Meredith.
}

\examples{
# uses the built-in volcano data set
require(grDevices) # for colours
require(graphics)
image(volcano, col=rampYOR(225), main="Sequential yellow-orange-red")
image(volcano, col=rampBYR(225), main="Divergent blue-yellow-red")
image(volcano, col=rampGBr(225), main="Divergent green-brown")
points(runif(500), runif(500), pch=16, cex=0.7) # add points
title(sub="Points are hard to see on the darkest colors", line=2)
image(volcano, col=rampGBr(225, range=2:10), main="Green-brown without darkest colors")
points(runif(500), runif(500), pch=16, cex=0.7)

# Try with a raster
require(raster)
r <- raster(system.file("external/test.grd", package="raster"))
plot(r, main="Default colors") # default is rev(terrain.colors(225))
plot(r, col=rampGBr(225), main="Divergent green-brown")
plot(r, col=rampGBr(225, bias=3), main="Divergent green-brown, bias=3")
plot(r, col=rampGBr(225, bias=0.5), main="Divergent green-brown, bias=0.5")
}
