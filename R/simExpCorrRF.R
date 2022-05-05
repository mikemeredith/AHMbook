# Define function for simulating spatially correlated random field
# AHM2 - 9.2

# Called by functions `simNmixSpatial` and `simOccSpatial`.

# Modified to use package 'fields' if 'RandomFields' is not available.

# In DESCRIPTION file:
# - add 'fields' to Imports
# - move RandomFields from Imports to Suggests

# In NAMESPACE file:
# - comment out or delete importFrom("RandomFields", "RFoptions", "RFsimulate", "RMexp")
# - add: importFrom("fields", "circulantEmbeddingSetup", "circulantEmbedding")



# ------------ Start of function definition ----------------
simExpCorrRF <- function(variance = 1, theta = 1, size = 50, seed = NA, show.plots = TRUE){
# Function creates Gaussian random field with negative
# exponential correlation and visualizes correlation and random field
#
# Function arguments:
# theta: parameter governing spatial correlation (=1/phi)
# ("large theta means long range of correlation")
# Note that RMexp is specified in terms of phi = 1/theta
# variance: variance of field, set at 1
# grid.size: Number of pixels in either direction
# show.plot: if TRUE, plots of the data will be displayed;
#  set to FALSE if you are running simulations or use inside of other fct's.

# Generate correlated random variables in a square
step <- 1
x <- seq(1, size, step)
y <- seq(1, size, step)
# grid <- as.matrix(expand.grid(x,y))
grid <- cbind(x = rep(x, each=size), y = y)

if(requireNamespace("RandomFields", quietly=TRUE)) {
  RandomFields::RFoptions(seed=seed)
  field <- RandomFields::RFsimulate(RandomFields::RMexp(var = variance, scale = theta),
      x=x, y=y, grid=TRUE)@data$variable1
  RandomFields::RFoptions(seed=NA)
} else {
 message("Using package 'fields' instead of 'RandomFields'; see help(simExpCorrRF).")
  if(!is.na(seed))
    set.seed(seed)  # Only for compatibility with RandomFields, better to set seed before calling simExpCommRF
  obj <- circulantEmbeddingSetup(grid=list(x=x, y=y), Covariance="Exponential", aRange=theta)
  tmp <- try(circulantEmbedding(obj), silent=TRUE)
  if(inherits(tmp, "try-error"))
    stop("Simulation of random field failed.\nTry with smaller values for 'size' or 'theta'.")
  field <- as.vector(tmp * sqrt(variance))
}

# Plots
# Correlation function
if(show.plots){
  oldpar <- par(mfrow = c(1,2), mar = c(5,5,4,2), "cex.main") ; on.exit(par(oldpar))
  tryPlot <- try( {
    dis <- seq(0.01, 20, by = 0.01)
    corr <- exp(-dis/theta)
    plot(dis, corr, type = "l", xlab = "Distance", ylab = "Correlation", ylim = c(0,1), col = "blue", lwd = 2)
    text(0.8*max(dis), 0.8, labels = paste("theta:", theta))

    # Random field
    # image(x, y, field,col=topo.colors(20), main = paste("Gaussian random field with \n negative exponential correlation (theta =", theta, ")"), cex.main = 1)
    par(mar = c(3,2,5,1))
    raster::plot(rasterFromXYZ(cbind(grid, field)), col=topo.colors(20),
    main = paste("Gaussian random field with \n negative exponential correlation (theta =", theta, ")"), cex.main = 1, legend=FALSE, box=FALSE)
    box()
  }, silent = TRUE)
  if(inherits(tryPlot, "try-error"))
    tryPlotError(tryPlot)
}

# Output
return(list(variance = variance, theta = theta, size = size, seed = seed,
  field = field,
  grid = grid))
} # ------------ End of function definition ----------------

