\name{e2dist}
\alias{e2dist}
\title{
Calculates the pair-wise distances between two sets of points
}
\description{
Calculates the euclidian distance between each of the points defined by the coordinates in 'x' and each of those in 'y'.
}
\usage{
e2dist(x, y)
}
\arguments{
  \item{x}{
a 2-column matrix or data frame with the x and y coordinates of a set of points.
}
  \item{y}{
a 2-column matrix or data frame with the x and y coordinates of a second set of points.
}
}
\value{
A \code{nrows(x)} x \code{nrows(y)} matrix with the pair-wise distances. 
}

\author{
J. Andrew Royle
}
\examples{

pts1 <- expand.grid(x = 1:5, y = 6:8)
pts2 <- cbind(x=runif(5, 1, 5), y=runif(5, 6, 8))
require(graphics)
plot(pts1)
points(pts2, pch=19, col='red')
e2dist(x=pts1, y=pts2)
# x and y can be the same:
e2dist(x=pts2, y=pts2)

}
