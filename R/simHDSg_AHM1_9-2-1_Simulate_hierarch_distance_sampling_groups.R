# Functions for the book Applied Hierarchical Modeling in Ecology (AHM)
# Marc Kery & Andy Royle, Academic Press, 2016.

# simHDSg - AHM1 section 9.2.1 p466

# Function to simulate data under HDS protocol with groups
#  (introduced in AHM1 Section 9.2.1)

simHDSg <- function(type = c("line", "point"), nsites = 100, lambda.group = 0.75, alpha0 = 0, alpha1 =
0.5, beta0 = 1, beta1 = 0.5, B = 4, discard0 = TRUE, show.plot=TRUE){
#
# Function simulates hierarchical distance sampling (HDS) data for groups under
# either a line (type = "line") or a point (type = "point") transect protocol
# and using a half-normal detection function (Buckland et al. 2001).
# Other function arguments:
# nsites: Number of sites (spatial replication)
# lambda.group: Poisson mean of group size
# alpha0, alpha1: intercept and slope of log-linear model relating sigma of
# half-normal detection function to group size
# beta0, beta1: intercept and slope of log-linear model relating the Poisson
# mean of the number of groups per unit area to habitat
# B: strip half width
#
# Checks and fixes for input data -----------------------------
nsites <- round(nsites[1])
stopifNegative(lambda.group, allowZero=FALSE)
stopifNegative(B, allowZero=FALSE)
# --------------------------------------------

type <- match.arg(type)

# Get covariates
habitat <- rnorm(nsites) # Simulated covariate
# Simulate abundance model for groups (Poisson GLM for N)
lambda <- exp(beta0 + beta1*habitat) # Density of groups per "square"
N <- rpois(nsites, lambda) # site-specific number of groups
N.true <- N # for point: inside of B
# Simulate observation model
data <- groupsize <- NULL
for(i in 1:nsites){
  if(N[i]==0){
    data <- rbind(data,c(i,NA,NA,NA,NA,NA)) # save site, y=1, u, v, d
    next
  }
  if(type=="line"){
    # Simulation of distances, uniformly, for each individual in the population
    d <- runif(N[i], 0, B)
    gs <- rpois(N[i],lambda.group) +1 # Observable group sizes >= 1
    groupsize <-c(groupsize,gs)
    sigma.vec <- exp(alpha0 + alpha1*(gs-1)) # Subtract 1 for interpretation
    # Detection probability for each group
    p <- exp(-d*d/(2*(sigma.vec^2)))
    # Determine if individuals are captured or not
    y <- rbinom(N[i], 1, p)
    u1 <- u2 <- rep(NA,N[i])
    # Subset to "captured" individuals only
    d <- d[y==1] ; u1 <- u1[y==1] ; u2 <- u2[y==1] ; gs <- gs[y==1] ; y <- y[y==1]
  }
  if(type=="point"){
    # Simulation of data on a circle of radius B (algorithm of Wallin)
    angle <- runif(N[i], 0, 2*pi)
    r2 <- runif(N[i], 0, 1)
    r <- B*sqrt(r2)
    u1 <- r*cos(angle) + B
    u2 <- r*sin(angle) + B

    d <- sqrt((u1 - B)^2 + (u2-B)^2)  ## d == r ! This block is all cruft
    N.true[i] <- sum(d<= B) # Population size inside of count circle, should be N[i] here.
    gs <- rpois(N[i], lambda.group) + 1
    groupsize <-c(groupsize,gs)
    sigma.vec <- exp(alpha0 + alpha1*(gs-1))
    # For counting individuals on a circle so we truncate p here ## cruft
    p <- ifelse(d<(B), 1, 0)*exp(-d*d/(2*(sigma.vec^2)))
    y <- rbinom(N[i], 1, p)
    # Subset to "captured" individuals only
    d <- d[y==1] ; u1 <- u1[y==1] ; u2 <- u2[y==1] ; gs <- gs[y==1] ; y <- y[y==1]
  }
  # Now compile things into a matrix and insert NA if no individuals were
  # captured at site i. Coordinates (u,v) are preserved.
  if(sum(y) > 0)  {
    data <- rbind(data,cbind(rep(i, sum(y)), y, u1, u2, d, gs))
  } else {
    data <- rbind(data,c(i,NA,NA,NA,NA,NA)) # make a row of missing data
  }
}
colnames(data)[1] <- "site"
# Subset to sites at which individuals were captured. You may or may not
# do this depending on how the model is formulated so be careful.
if(discard0)
  data <- data[!is.na(data[,2]),]
# Visualization
if(show.plot) {
  if(type=="line"){ # For line transect
    op <- par(mfrow = c(1, 3)) ; on.exit(par(op))
    tryPlot <- try( {
      hist(data[,"d"], col = "lightblue", breaks = 20, main =
      "Frequency of distances to groups", xlab = "Distance")
      ttt <- table(data[,1])
      n <- rep(0, nsites)
      n[as.numeric(rownames(ttt))] <- ttt
      plot(habitat, n, main = "Observed group counts (n) vs. habitat", frame = FALSE)
      plot(table(data[,"gs"]), main = "Observed group sizes", ylab = "Frequency", frame = FALSE)
    }, silent = TRUE)
    if(inherits(tryPlot, "try-error"))
      tryPlotError(tryPlot)
  }
  if(type=="point"){ # For point transect
    op <- par(mfrow = c(2,2)) ; on.exit(par(op))
    tryPlot <- try( {
      plot(data[,"u1"], data[,"u2"], pch = 16, main =
      "Located groups in point transects", xlim = c(0, 2*B),
      ylim = c(0, 2*B), col = data[,1], asp = 1)
      points(B, B, pch = "+", cex = 3)
      plotrix::draw.circle(B, B, B)
      hist(data[,"d"], col = "lightblue", breaks = 20, main =
      "Frequency of distances to groups", xlab = "Distance")
      ttt <- table(data[,1])
      n <- rep(0, nsites)
      n[as.numeric(rownames(ttt))] <- ttt
      plot(habitat, n, main = "Observed group counts (n) vs. habitat", frame = FALSE)
      plot(table(data[,"gs"]), main = "Observed group sizes", ylab = "Frequency", frame = FALSE)
    }, silent = TRUE)
    if(inherits(tryPlot, "try-error"))
      tryPlotError(tryPlot)
  }
}

# Output
list(type = type, nsites = nsites, lambda.group = lambda.group, alpha0 = alpha0,
  alpha1 = alpha1, beta0 = beta0, beta1 = beta1, B = B, data=data, habitat=habitat,
  N = N, N.true = N.true, groupsize=groupsize)
}


