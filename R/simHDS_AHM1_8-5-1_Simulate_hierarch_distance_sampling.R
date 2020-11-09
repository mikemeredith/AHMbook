# Functions for the book Applied Hierarchical Modeling in Ecology (AHM)
# Marc Kery & Andy Royle, Academic Press, 2016.

# simHDS  - AHM1 section 8.5.1 p444

# Function to simulate data under hierarchical distance sampling protocol (line or point)
#   (introduced in AHM1 Section 8.5.1)
simHDS  <- function(type=c("line", "point"), nsites = 100, mean.lambda = 2,
   beta.lam  = 1, mean.sigma = 1, beta.sig = -0.5, B = 3, discard0=TRUE, show.plot=TRUE){
#
# Function simulates hierarchical distance sampling (HDS) data under
#   either a line (type = "line") or a point (type = "point")
#   transect protocol.
#   Function arguments:
#     nsites: Number of sites (spatial replication)
#     alpha.lam (= log(mean.lambda)), beta.lam: intercept and
#        slope of log-linear regression of expected lambda
#        on a habitat covariate
#     alpha.sig (= log(mean.sigma)), beta.sig: intercept and
#        slope of log-linear regression of scale parameter of
#        half-normal detection function on wind speed
#     B: strip half width

# Checks and fixes for input data -----------------------------
nsites <- round(nsites[1])
stopifNegative(mean.lambda, allowZero=FALSE)
stopifNegative(mean.sigma, allowZero=FALSE)
stopifNegative(B, allowZero=FALSE)
# --------------------------------------------

type <- match.arg(type)

# Get covariates
habitat <- rnorm(nsites)                    # habitat covariate
wind <- runif(nsites, -2, 2)                # wind covariate

# Simulate abundance model (Poisson GLM for N)
lambda <- exp(log(mean.lambda) + beta.lam*habitat) # density per "square"
N <- rpois(nsites, lambda)                  # site-specific abundances
N.true <- N                                 # for point: inside B

# Detection probability model (site specific)
sigma <- exp(log(mean.sigma) + beta.sig*wind)

# Simulate observation model
data <- NULL

for(i in 1:nsites){
  if(N[i]==0){
    data <- rbind(data, c(i,NA,NA,NA,NA)) # save site, y=1, u, v, d
  next
  }
  if(type=="line"){
    # Simulation of distances, uniformly, for each ind. in pop.
    # note it piles up all N[i] guys on one side of the transect
    d <- runif(N[i], 0, B)
    p <- exp(-d *d / (2 * (sigma[i]^2)))
    # Determine if individuals are captured or not
    y <- rbinom(N[i], 1, p)
    u <- v <- rep(NA, N[i])   # coordinates (u,v)
    # Subset to "captured" individuals only
    d <- d[y==1]
    u <- u[y==1]
    v <- v[y==1]
    y <- y[y==1]
  }

  if(type=="point"){
    # Simulation data on a square
    u <- runif(N[i], 0, 2*B)
    v <- runif(N[i], 0, 2*B)
    d <- sqrt((u-B)^2 + (v-B)^2)
    N.true[i] <- sum(d<= B)    # Population size inside of count circle

    # Can only count indidividuals in the circle, so set to zero p
    # of individuals in the corners (thereby truncating them)
    p <- exp(-d *d / (2 * (sigma[i]^2)))  # Detection probabiilty ..
    pp <- ifelse(d <= B, 1, 0) * p    # ... times "inside or outside"
    y <- rbinom(N[i], 1, pp)  # Det./non-detection of each individual

    # Subset to "captured" individuals only
    u <- u[y==1]
    v <- v[y==1]
    d <- d[y==1]
    y <- y[y==1]
  }

  # Compile things into a matrix and insert NA if no individuals were
  # captured at site i. Coordinates (u,v) are not used here.
  if(sum(y) > 0)
    data <- rbind(data, cbind(rep(i, sum(y)), y, u, v, d))
  else
    data <- rbind(data, c(i,NA,NA,NA,NA)) # make a row of missing data
}
colnames(data) <- c("site", "y", "u", "v", "d") # name 1st col "site"

# Subset to sites at which individuals were captured. You may or may not
#  want to do this depending on how the model is formulated so be careful.
if(discard0)
  data <- data[!is.na(data[,2]),]

# Visualisation
if(show.plot) {
  if(type=="line"){       # For line transect
    op <- par(mfrow = c(1, 3)) ; on.exit(par(op))
    tryPlot <- try( {
      hist(data[,"d"], col = "lightblue", breaks = 20, main =
        "Frequency of distances", xlab = "Distance")
      ttt <- table(data[,1])
      n <- rep(0, nsites)
      n[as.numeric(rownames(ttt))] <- ttt
      plot(habitat, n, main = "Observed counts (n) vs. habitat")
      plot(wind, n, main = "Observed counts (n) vs. wind speed")
    }, silent = TRUE)
    if(inherits(tryPlot, "try-error"))
      tryPlotError(tryPlot)
  }

  if(type=="point"){       # For point transect
    op <- par(mfrow = c(2,2)) ; on.exit(par(op))
    tryPlot <- try( {
      plot(data[,"u"], data[,"v"], pch = 16, main =
        "Located individuals in point transects", xlim = c(0, 2*B),
        ylim = c(0, 2*B), col = data[,1], asp = 1)
      points(B, B, pch = "+", cex = 3, col = "black")
      plotrix::draw.circle(B, B, B)
      hist(data[,"d"], col = "lightblue", breaks = 20, main =
        "Frequency of distances", xlab = "Distance")
      ttt <- table(data[,1])
      n <- rep(0, nsites)
      n[as.numeric(rownames(ttt))] <- ttt
      plot(habitat, n, main = "Observed counts (n) vs. habitat")
      plot(wind, n, main = "Observed counts (n) vs. wind speed")
    }, silent = TRUE)
    if(inherits(tryPlot, "try-error"))
      tryPlotError(tryPlot)
  }
}

# Output
list(type = type, nsites = nsites, mean.lambda = mean.lambda, beta.lam = beta.lam,
    mean.sigma = mean.sigma, beta.sig = beta.sig, B = B, data=data, habitat=habitat,
    wind=wind, N = N, N.true = N.true )
}

