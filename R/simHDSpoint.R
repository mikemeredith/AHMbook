# Functions for the book Applied Hierarchical Modeling in Ecology (AHM)
# Marc Kery & Andy Royle, Academic Press, 2016.

# simHDSpoint

# This is a re-factored and more efficient version of the 'point' type for 'simHDS'.
# Note that some of the arguments have changed.

# Function to simulate data under hierarchical distance sampling protocol point transects

simHDSpoint  <- function(nsites = 1000, mean.density = 1,
   beta.density  = 1, mean.sigma = 20, beta.sigma = -5, B = 60, discard0=FALSE, show.plots=TRUE){
#
# Function simulates hierarchical distance sampling (HDS) data under
#   a point transect protocol.

#  Function arguments:
#     nsites: Number of sites (spatial replication)
#     mean.density: expected DENSITY, number of individuals per HECTARE
#     beta.density: coefficient of log(expected density) on habitat covariate
#     mean.sigma: scale parameter sigma in the half-normal detection function in METERS
#     beta.sigma: slope of log-linear regression of scale parameter on wind speed
#     B: maximum distance from the point, in METERS.
#     discard0: if TRUE, sites with no detections will be removed from the output.
#     show.plots: if TRUE, plots of the output are displayed

# Checks and fixes for input data -----------------------------
nsites <- round(nsites[1])
stopifNegative(mean.density, allowZero=FALSE)
stopifNegative(mean.sigma, allowZero=FALSE)
stopifNegative(B, allowZero=FALSE)
# --------------------------------------------

# Get covariates
habitat <- rnorm(nsites)                    # habitat covariate
wind <- runif(nsites, -2, 2)                # wind covariate

# Simulate abundance model (Poisson GLM for N)
lambda <- exp(log(mean.density) + beta.density*habitat) * base::pi * B^2 / 1e4 # expected number in circle
N <- rpois(nsites, lambda)                  # site-specific abundances

# Detection probability model (site specific)
sigma <- exp(log(mean.sigma) + beta.sigma*wind)

# Simulate observation model
dataList <- Nlist <- vector(mode = "list", length = nsites)
counts <- numeric(nsites)

for(i in 1:nsites){
  if(N[i]==0){
    dataList[[i]] <- c(site=i, d=NA)
    next
  }
  # Simulation of distance from point given uniform distribution on the circle (algorithm of Wallin)
  d <- sqrt(runif(N[i], 0, 1)) * B
  Nlist[[i]] <- d

  # Detection process, half-normal detection function
  p <- exp(-d^2 / (2 * (sigma[i]^2)))  # Detection probability ..
  det <- rbinom(N[i], 1, p)              # Det./non-detection of each individual
  counts[i] <- sum(det)
  # Subset to "captured" individuals only
  d <- d[det==1]
  if(length(d) == 0)
    d <- NA

  # Plug the d's into the list together with site IDs.
  dataList[[i]] <- cbind(site=i, d=d)
}

# convert dataList to a single matrix
data <- do.call(rbind, dataList)

# Subset to sites at which individuals were captured. You rarely
#  want to do this depending on how the model is formulated so be careful.
if(discard0)
  data <- data[!is.na(data[,2]),]

# Visualisation
if(show.plots) {
  op <- par(mfrow = c(2,2)) ; on.exit(par(op))
  tryPlot <- try( {

    angle <- runif(nrow(data), 0, 2*base::pi)
    x <- data[, 'd'] * cos(angle)
    y <- data[, 'd'] * sin(angle)
    plot(x, y, pch = 16, main =
        "Locations of detected individuals", xlim = c(-B, B),
        ylim = c(-B, B), col = data[,1], asp = 1)
    points(0, 0, pch = "+", cex = 3, col = "black")
    plotrix::draw.circle(0, 0, B)

    tmp <- hist(unlist(Nlist), breaks = 20, xlim=c(0, B), xlab = "Distance",
        main = "Frequency of distances\nblue: detected, gray: undetected")
    hist(data[,"d"], breaks = tmp$breaks,
        col = "lightblue", add=TRUE)

    plot(habitat, counts, main = "Observed counts (n) vs. habitat")

    plot(wind, counts, main = "Observed counts (n) vs. wind speed")
  }, silent = TRUE)
  if(inherits(tryPlot, "try-error"))
    tryPlotError(tryPlot)
  }


# Output
list( # arguments input
    nsites = nsites, mean.density = mean.density, beta.density = beta.density,
    mean.sigma = mean.sigma, beta.sigma = beta.sigma, B = B,
    # generated values
    data = data, counts = counts, habitat = habitat, wind = wind, N = N )
}

