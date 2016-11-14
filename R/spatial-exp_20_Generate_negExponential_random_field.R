# Functions for the book Applied Hierarchical Modeling in Ecology (AHM)
# Marc Kéry & Andy Royle, Academic Press, 2016.

# spatial.exp - section 20.??


# Function to generate random field with negative exponential correlation

spatial.exp <- function(variance = 1, theta = 1, size = 50, show.plot = TRUE){
# Function creates Gaussian random field with negative
# exponential correlation and visualizes correlation and random field
#
# Function arguments:
# theta: parameter governing spatial correlation (=1/phi)
# ("large theta means high correlation")
# Note that RMexp is specified in terms of phi = 1/theta
# variance: variance of field, set at 1
# grid.size: Number of pixels in either direction
# show.plot: if TRUE, plots of the data will be displayed;
#  set to FALSE if you are running simulations or use inside of other fct's.

# library(raster)
# library(RandomFields)

# Generate correlated random variables in a square
RFoptions(seed=NA)
step <- 1
size<-size
x <- seq(1, size, step)
y <- seq(1, size, step)
grid <- as.matrix(expand.grid(x,y))
field <- matrix(RFsimulate(RMexp(var = variance, scale = theta), x=x, y=y, grid=TRUE)@data$variable1, ncol = size)

# Plots
# Correlation function
if(show.plot){
par(mfrow = c(1,2), mar = c(5,5,4,2))
dis <- seq(0.01, 20, by = 0.01)
corr <- exp(-dis/theta)
plot(dis, corr, type = "l", xlab = "Distance", ylab = "Correlation", ylim = c(0,1), col = "blue", lwd = 2)
text(0.8*max(dis), 0.8, labels = paste("theta:", theta))

# Random field
image(x, y, field,col=topo.colors(20), main = paste("Gaussian random field with \n negative exponential correlation (theta =", theta, ")"), cex.main = 1)
}

# Output
return(list(variance = variance, theta = theta, field = field, grid = grid))
}

