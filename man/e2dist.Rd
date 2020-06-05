\encoding{UTF-8}
\name{e2dist}
\alias{e2dist}
\title{
Calculates the pair-wise distances between two sets of points
}
\description{
Calculates the Euclidian distance between each of the points defined by the coordinates in 'x' and each of those in 'y'.
}
\usage{
e2dist(x, y = NULL)
}
\arguments{
  \item{x}{
a 2-column matrix or data frame with the x and y coordinates of a set of points, or a length 2 vector with the coordinates of a single point.
}
  \item{y}{
an optional 2-column matrix or data frame, or a length 2 vector, with the x and y coordinates of a second set of points; if NULL, y is taken to be the same as x.
}
}
\value{
A \code{nrows(x)} x \code{nrows(y)} matrix with the pair-wise distances.
}

\author{
Andy Royle
}
\examples{

pts1 <- expand.grid(x = 1:5, y = 6:8)
pts2 <- cbind(x=runif(5, 1, 5), y=runif(5, 6, 8))
require(graphics)
plot(pts1)
points(pts2, pch=19, col='red')
e2dist(x=pts1, y=pts2)
# with a vector argument:
vec <- colMeans(pts2)
e2dist(vec, pts2)
# with y = NULL:
e2dist(x=pts2)
}
