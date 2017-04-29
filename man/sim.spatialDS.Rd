\name{sim.spatialDS}
\alias{sim.spatialDS}
\encoding{UTF-8}

\title{
Simulates data for a basic spatial distance sampling model
}
\description{
Generates data with the following steps:

1. Simulate a spatially correlated habitat covariate (x) over a grid of pixels covering a square.

2. Distribute the population of N individuals over the square with probability of location in a pixel related to the covariate.

3. Decide which individuals are detected using a distance sampling model with an observer at the centre of the square, with either a half normal or a hazard rate detection function. (Note that individuals outside the circle of radius \code{B} can be detected.)

4. If \code{keep.all = FALSE}, return the locations of only the individuals detected.
}
\usage{
sim.spatialDS(N = 1000, beta = 1, sigma = 1, keep.all = FALSE, B = 3,
  model=c("logit", "halfnorm", "hazard"), lambda = B/3, show.plot=TRUE)
}
\arguments{
  \item{N}{
total population size in the square
}
  \item{beta}{
coefficient of SOMETHING on spatial covariate x
}
  \item{sigma}{
scale parameter of detection function
}
  \item{keep.all}{
if TRUE, the data for all individuals are returned; if FALSE, only for individuals detected.
}
  \item{B}{
circle radius; the data are simulated over a square of side 2 * B. This is usually set so that the probability of detection of individuals outside the circle is negligable, eg, \code{B = 3*sigma}. %, but individuals outside the circle of radius B are not detected. % In this case they ARE detected.
}
  \item{model}{
The detection function used, can be "logit", "halfnorm" or "hazard": see Details.
}
  \item{lambda}{
The scale parameter for the spatially autocorrelated Habitat covariate.
}
  \item{show.plot}{
choose whether to show plots or not. Set to FALSE when using function in simulations.
}
}
\details{
The "logit" detection function is \code{2*plogis(-d^2/(2*sigma^2))}, which corresponds to the detection model implemented in \code{unmarked::pcount.spHDS}.

The hazard rate detection model used here is \code{p <- 1-exp(-exp(-d^2/(2*sigma^2)))}. This has no shape parameter, and probability of detection when \code{d = 0} is \code{p(0) = 1-exp(-exp(0)) = 0.632}. This is chosen so that the complementary log-log of p has only a quadratic term of distance: \code{-(1/(2*sigma^2)) * d^2}. See Kéry & Royle (2016) p.540 for details.
}
\value{
A list with the values of the input arguments and the following additional elements:
  \item{u1 }{x coordinate of each animal}
  \item{u2 }{y coordinate of each animal}
  \item{d }{distance of each animal from the centre of the circle}
  \item{pixel.id }{the pixel in which each animal is located, the row number in \code{grid}}
  \item{y }{indicator of detection of each animal, a vector of length N}
  \item{N.real }{the number of animals inside the circle of radius B}
  \item{Habitat }{Value of the spatially  correlated habitat covariate, a 900 x 1 matrix}
  \item{grid }{Coordinates of the centres of each pixel, a dataframe with 900 rows and 2 columns}
If \code{keep.all = FALSE} (the default), only the animals detected are included in \code{u1, u2, d, pixel.id}.
}
\references{
Kéry & Royle (2016) \emph{Applied Hierarachical Modeling in Ecology} Section 9.8.3

}
\author{
Kéry & Royle
}

\examples{
# Generate data with the default arguments and look at the structure:
tmp <- sim.spatialDS()
str(tmp)

# Generate date with model = "logit" and analyse the data with unmarked::pcount.spatialHDS
set.seed(1234)
tmp <- sim.spatialDS(model="logit")
# Plot shows a large area of good habitat west of the observer with many animals detected
str(tmp)  # 272 animals detected out of 850 inside the circle (N.real)

# Get the count of animals detected in each pixel
pixel.count <- tabulate(tmp$pixel.id, nbins=nrow(tmp$grid))
# Centre the Habitat covariate
Habitat <- tmp$Habitat - mean(tmp$Habitat)
# Create a detection covariate: distance between observer and pixel center
dist <- with(tmp, sqrt((grid[,1]-B)^2 + (grid[,2]-B)^2))
# Construct an unmarkedFrame
umf <- unmarkedFramePCount(y=cbind(pixel.count),
   siteCovs=data.frame(dist=dist, Habitat=Habitat))
summary(umf)

# Fit some models
(fm0 <- pcount.spHDS(~ -1 + I(dist^2) ~ 1, umf, K=20))
(fm1 <- pcount.spHDS(~ -1 + I(dist^2) ~ Habitat, umf, K=20))
# The model with Habitat has much lower AIC
# Get an estimate of the total population in the square (true is N = 1000)
sum(predict(fm1, type='state')[, 1])
}
