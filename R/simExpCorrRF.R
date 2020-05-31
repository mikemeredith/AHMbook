# Define function for simulating spatially correlated random field

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

# library(raster)
# library(RandomFields)

# Generate correlated random variables in a square
step <- 1
x <- seq(1, size, step)
y <- seq(1, size, step)
# grid <- as.matrix(expand.grid(x,y))
grid <- cbind(x = rep(x, each=size), y = y)

RFoptions(seed=seed)
# field <- matrix(RFsimulate(RMexp(var = variance, scale = theta), x=x, y=y, grid=TRUE)@data$variable1, ncol = size)
field <- RFsimulate(RMexp(var = variance, scale = theta), x=x, y=y, grid=TRUE)@data$variable1
RFoptions(seed=NA)

# Plots
# Correlation function
if(show.plots){
  oldpar <- par(mfrow = c(1,2), mar = c(5,5,4,2), "cex.main") ; on.exit(par(oldpar))
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
}

# Output
return(list(variance = variance, theta = theta, size = size, seed = seed,
  field = field,
  grid = grid))
} # ------------ End of function definition ----------------


