# Functions for the book Applied Hierarchical Modeling in Ecology (AHM)
# Marc Kery & Andy Royle, Academic Press, 2016.

# sim.ldata  - AHM1 section 8.2.3 p402



# Function to simulate non-hierarchical line transect data
#   (introduced in AHM1 Section 8.2.3)
sim.ldata <- function(N = 200, sigma = 30, show.plot = TRUE){
# Function to simulate line transect data under CDS.
# Function arguments:
#    N: number of individuals along transect with distance u(-100, 100)
#    sigma: scale parameter of half-normal detection function
# Function subjects N individuals to sampling, and then retains the value
# of x=distance only for individuals that are captured
if(FALSE) x <- NULL # fix issues with curve

# Checks and fixes for input data -----------------------------
N <- round(N[1])
stopifNegative(sigma, allowZero=FALSE)
# --------------------------------------------

xall <- runif(N, -100,100) # Distances of all N individuals
g <- function(x, sig) exp(-x^2/(2*sig^2))
p <- g(xall, sig=sigma) # detection probability
y <- rbinom(N, 1, p) # some inds. are detected and their distance measured
x <- xall[y==1]      # this has direction (right or left transect side)
x <- abs(x)          # now it doesn't have direction
if(show.plot) {
  op <- par(mfrow = c(1,2)) ; on.exit(par(op))
  # Plot the detection function
  tryPlot <- try( {
    curve(exp(-x^2/(2*sigma^2)), 0, 100, xlab="Distance (x)", ylab="Detection prob.",
        lwd = 2, main = "Detection function", ylim = c(0,1))
    text(80, 0.9, paste("sigma:", sigma))
    hist(abs(xall), nclass=10, xlab = "Distance (x)", col = "grey", 
        main = "True (grey) \nand observed distances (blue)")
    hist(x, col = "blue", add = TRUE)
  }, silent = TRUE)
  if(inherits(tryPlot, "try-error"))
    tryPlotError(tryPlot)
}
return(list(N = N, sigma = sigma, xall = xall, x = x))
}


