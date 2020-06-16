# Functions for the book Applied Hierarchical Modeling in Ecology (AHM)
# Marc Kery & Andy Royle, Academic Press, 2016.

# sim.pdata  - AHM1 section 8.2.5.1 p410

# Function to simulate non-hierarchical point transect (= point count) data
#   (introduced in AHM1 Section 8.2.5.1)

sim.pdata <- function(N=1000, sigma=1, B=3, keep.all=FALSE, show.plot=TRUE) {
# Function simulates coordinates of individuals on a square
# Square is [0,2*B] x[0,2*B], with a count location on the center
# point (B,B)
# Function arguments:
#    N: total population size in the square
#    sigma: scale of half-normal detection function
#    B: circle radias
#    keep.all: return the data for y = 0 individuals or not
if(FALSE) x <- NULL # Kludge to keep R CMD check happy with curve

# Checks and fixes for input data -----------------------------
N <- round(N[1])
stopifNegative(sigma, allowZero=FALSE)
stopifNegative(B, allowZero=FALSE)
# --------------------------------------------

# Simulate and plot simulated data
u1 <-runif(N, 0, 2*B)           # (u1,u2) coordinates of N individuals
u2 <- runif(N, 0, 2*B)
d <- sqrt((u1 - B)^2 + (u2 - B)^2) # distance to center point of square
N.real <- sum(d<= B)           # Population size inside of count circle

# Can only count indidividuals in the circle, so set to zero detection probability of individuals in the corners (thereby truncating them):
p <- ifelse(d < B, 1, 0) * exp(-d*d/(2*(sigma^2)))
# Now we decide whether each individual is detected or not
y <- rbinom(N, 1, p)

if(show.plot) {
  op <- par(mfrow = c(1,2)) ; on.exit(par(op))
  # Plot the detection function
  tryPlot <- try( {
    curve(exp(-x^2/(2*sigma^2)), 0, B, xlab="Distance (x)", ylab="Detection prob.",
      lwd = 2, main = "Detection function", ylim = c(0,1))
    text(0.8*B, 0.9, paste("sigma:", sigma))
    plot(u1, u2, asp = 1, pch = 1, main = "Point transect")
    points(u1[d <= B], u2[d <= B], pch = 16, col = "black")
    points(u1[y==1], u2[y==1], pch = 16, col = "blue")
    points(B, B, pch = "+", cex = 3, col = "red")
    plotrix::draw.circle(B, B, B)
  }, silent = TRUE)
  if(inherits(tryPlot, "try-error"))
    tryPlotError(tryPlot)
}

# Put all of the data in a matrix:
#      (note we don't care about y, u, or v normally)

if(!keep.all){
   u1 <- u1[y==1]
   u2 <- u2[y==1]
   d <- d[y==1]
}
return(list(N=N, sigma=sigma, B=B, u1=u1, u2=u2, d=d, y=y, N.real=N.real))
}
