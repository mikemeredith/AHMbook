\name{simPPe}
\alias{simPPe}
\encoding{UTF-8}

\title{
Simulate a spatial point pattern in a heterogeneous landscape
}
\description{
The function simulates a spatial point pattern in a heterogeneous landscape simulated on a square landscape. The study area ('core') is simulated inside the larger landscape that includes a buffer. The size of the core is defined by the lscape.size minus twice the buffer.

There is one habitat covariate X that affects the intensity of the points. X is spatially structured with negative-exponential spatial autocorrelation; the parameters of the field can be chosen to create large 'islands' of similar values or no 'islands' at all, in which case the field is spatially unstructured.

The intensity of STATIC points (e.g. animal activity centers) may be inhomogeneous and affected by the coefficient beta, which is the log-linear effect of X.
}
\usage{
simPPe(lscape.size = 150, buffer.width = 25, variance.X = 1, theta.X = 10,
  M = 250, beta = 1, quads.along.side = 6)
}

\arguments{
  \item{lscape.size}{
size (width = height) of the square landscape, including core and buffer.
}
  \item{buffer.width}{
width of buffer around core study area.
}
  \item{variance.X}{
variance of Gaussian random field (covariate X).
}
  \item{theta.X}{
scale parameter of correlation in field (must be >0).
}
  \item{M}{
expected number of activity centers in core area.
}
  \item{beta}{
coefficient of the habitat covariate.
}
  \item{quads.along.side}{
number of quadrats along the side of the core area; the total number of quadrats will be quads.along.side^2, thus indirectly defining the quadrat area.
}
}

\value{
A list with the values of the input arguments and the following additional elements:
 \item{core }{range of x and y coordinates in the 'core'}
 \item{M2 }{number of ACs in the total landscape, including the buffer}
 \item{grid }{coordinates of the centre of each pixel}
 \item{pixel.size }{length of side of each pixel}
 \item{size.core }{the width=height of the core area}
 \item{prop.core }{the proportion of the landscape inside the core}
 \item{X }{matrix of covariate values for each pixel}
 \item{probs }{matrix of probabilities of an AC being inside each pixel (sums to 1)}
 \item{pixel.id }{the ID of the pixel for each AC}
 \item{u }{2-column matrix, coordinate of each AC}
 \item{nsite }{number of quadrats}
 \item{quad.size }{width = height of each quadrat}
 \item{breaks }{boundaries of the quadrats}
 \item{mid.pt }{mid-points of the quadrats}
 \item{lambda_pp }{intensity of point pattern (ACs per unit area)}
 \item{Nac }{site-specific abundance of ACs}
 \item{zac }{site-specific occurrence (0/1) of ACs}
 \item{E_N }{average realized abundance per quadrat}
 \item{E_z }{average realized occupancy per quadrat}
}
\references{
Kéry, M. & Royle, J.A. (2020) \emph{Applied Hierarchical Modeling in Ecology} Chapter 10
}
\author{
Marc Kéry & Andy Royle.
}

\examples{
# Nice plot (produces the really nice Fig. 10-2 in the book)
set.seed(117)
str(dat <- simPPe(lscape.size = 200, buffer.width = 25, variance.X = 1,
  theta.X = 70, M = 200, beta = 1, quads.along.side = 6))

\donttest{
# Smaller study area, fewer individuals (M)
str(dat <- simPPe(lscape.size = 24, buffer.width = 2, variance.X = 1,
  theta.X = 10, M = 50, beta = 1, quads.along.side = 6))

# Stronger habitat heterogeneity (variance.X): more aggregation
str(dat <- simPPe(lscape.size = 24, buffer.width = 2, variance.X = 10,
  theta.X = 10, M = 50, beta = 1, quads.along.side = 6))

# Longer habitat gradient (theta.X)
str(dat <- simPPe(lscape.size = 24, buffer.width = 2, variance.X = 1,
  theta.X = 100, M = 250, beta = 1, quads.along.side = 6))

# No habitat variability (variance.X): homogeneous point process
str(dat <- simPPe(lscape.size = 24, buffer.width = 2, variance.X = 0,
  theta.X = 10, M = 100, beta = 1, quads.along.side = 6))
# No habitat preference (beta): homogeneous point process
str(dat <- simPPe(lscape.size = 24, buffer.width = 2, variance.X = 1,
  theta.X = 10, M = 100, beta = 0, quads.along.side = 6))
# Habitat heterogeneity at very small scale (theta.X)
str(dat <- simPPe(lscape.size = 1000, buffer.width = 20, variance.X = 1,
  theta.X = 0.001, M = 250, beta = 1, quads.along.side = 6))

str(simPPe(M = 1))      # Often produces no point at all
str(simPPe(M = 10))
str(simPPe(M = 100))
str(simPPe(M = 1000))

str(simPPe(M = 20, quads.along.side = 100)) # Lots of small sites
str(simPPe(M = 20, quads.along.side = 10))
str(simPPe(M = 20, quads.along.side = 5))
str(simPPe(M = 20, quads.along.side = 1))   # whole study area is 1 site
}
}
