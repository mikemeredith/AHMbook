# Functions for the book Applied Hierarchical Modeling in Ecology (AHM)
# Marc Kery & Andy Royle, Academic Press, 2016.

# sim.spatialDS - AHM1 section 9.8.3 p534

# Function generates data under spatial hierarchical distance sampling model
#   (introduced in AHM1 Section 9.8.3)

# Changes April/May 2017 (v. 0.1.3.9000/9001)
#   Added lambda
#   Added logit detection function, removed hazard
#   Added useHabitat option to simulate new detection data with old Habitat values

sim.spatialDS <- function(N=1000, beta = 1, sigma=1, keep.all=FALSE,
  B=3, model=c("logit", "halfnorm"), lambda = B/3, useHabitat, show.plot=TRUE){
# Function simulates coordinates of individuals on a square
# Square is [0,2B] x [0,2B], with a count location on the point (B, B)
#   N: total population size in the square
#   beta: coefficient of SOEMTHING on spatial covariate x
#   sigma: scale of half-normal detection function
#   keep.all: return the data for all individuals, whether detected or not

# Checks and fixes for input data -----------------------------
N <- round(N[1])
stopifNegative(sigma, allowZero=FALSE)
stopifNegative(B, allowZero=FALSE)
# --------------------------------------------

model <- match.arg(model)

# Create coordinates for 30 x 30 grid
delta <- (2*B-0)/30                # '2D bin width'
grx <- seq(delta/2, 2*B - delta/2, delta) # mid-point coordinates
gr <- expand.grid(grx, grx, KEEP.OUT.ATTRS = FALSE) # Create grid coordinates

if(missing(useHabitat)) {
  # Create spatially correlated covariate x
  V <- exp(-e2dist(gr,gr)/lambda)
  x <- t(chol(V))%*%rnorm(900)
} else {
  x <- useHabitat$Habitat
  # Check input
  if(is.null(x) ||
      is.null(dim(x)) ||
      dim(x)[2] != 1 ||
      dim(x)[1] != 900)
    stop("useHabitat is not valid output from sim.spatialDS.")
}

# Simulate point locations as function of habitat covariate x
probs <- exp(beta*x)/sum(exp(beta*x)) # probability of point in pixel (sum = 1)
pixel.id <- sample(1:900, N, replace=TRUE, prob=probs)
# could simulate randomly within the pixel but it won't matter so place centrally
u1 <- gr[pixel.id,1]
u2 <- gr[pixel.id,2]

d <- sqrt((u1 - B)^2 + (u2-B)^2)   # distance to center point of square
#plot(u1, u2, pch = 1, main = "Point transect")
N.real <- sum(d <= B)               # Population size inside of count circle

# Can only count individuals in the circle, so set to zero detection probability of individuals in the corners (thereby truncating them)
# p <- ifelse(d < B, 1, 0) * exp(-d*d/(2*(sigma^2)))
# We do away with the circle constraint here.
if(model=="halfnorm")
   p <- exp(-d*d/(2*sigma*sigma))
if(model=="logit")
   p<- 2*plogis(  -d*d/(2*sigma*sigma)   )
# Now we decide whether each individual is detected or not
y <- rbinom(N, 1, p)                                           # detected or not

if(show.plot) {
  op <- par(mar=c(3,3,3,6)) ; on.exit(par(op))
  tryPlot <- try( {
    image(rasterFromXYZ(cbind(as.matrix(gr),x)), col=topo.colors(10), asp=1, bty='n') # need to convert gr to a matrix
    rect(0, 0, 2*B, 2*B)  # draw box around the image
    # draw.circle(B, B, B)
    points(B, B, pch="+", cex=3)
    image_scale(x, col=topo.colors(10))
    # points(u1, u2, pch=20, col='black', cex = 0.8)  # plot points
    title("Extremely cool figure")         # express your appreciation of all this
    # points(u1[d <= B], u2[d <= B], pch = 16, col = "black", cex = 1) # in circle but not detected
    points(u1, u2, pch = 16, col = c("black", "red")[y+1]) # in circle but not detected
    # points(u1[y==1], u2[y==1], pch = 16, col = "red", cex = 1)     # detected
  }, silent = TRUE)
  if(inherits(tryPlot, "try-error"))
    tryPlotError(tryPlot)
}
# Remove data for individuals not detected
if(!keep.all){
   u1 <- u1[y==1]
   u2 <- u2[y==1]
   d <- d[y==1]
   pixel.id <- pixel.id[y==1]
}
# Output
return(list(model=model, N=N, beta=beta, B=B, u1=u1, u2=u2, d=d, pixel.id=pixel.id,
  y=y, N.real=N.real, Habitat=x, grid=gr))
}

