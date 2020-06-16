# Functions for the book Applied Hierarchical Modeling in Ecology (AHM)
# Marc Kery & Andy Royle, Academic Press, 2016.

# From Marc, 13 June 2019
# "Build up the BerneseOberland data set and the simNmixSp functions.docx"

# Originally simNmixSpatial.

# ------  Define function simNmixSpatial --------
simNmixSpatial <- function(nsurveys = 3, mean.lambda = exp(2), beta = c(2, -2),
  mean.p = 0.5, alpha = c(-1, -1), sample.size = 500, variance.RF = 1, theta.RF = 10,
  seeds = c(10, 100), truncN = 6, show.plots = TRUE, verbose = TRUE) {

# Simulates replicated counts under a spatial, static binomial N-mixture model for a semi-realistic landscape in a square of 50x50 km in the BernesOberland around Interlaken, Switzerland.
# Unit of the data simulation is a 1km2 quadrat, hence, there are 2500 units.
# For abundance, the function allows you to specify a quadratic effect of elevation, the data for which are contained in the data set BerneseOberland, which is part of the AHMbook package and is a subset of the data set 'Switzerland' in R package unmarked.
# Then, a Gaussian spatial random field (s) with negative exponential correlation function is simulated using the AHMbook function simExpCorrRF. For that field, you can set the variance and the range scale parameter theta (see helptext for that function for more details). Basically, the larger the value of theta.RF, the bigger are the 'islands' simulated in the random field.
# The abundance in each quadrat i is built up via the following linear predictor:
#  lambda[i] <- exp(beta0 + beta1 * elev[i] + beta2 * elev[i]^2 + s[i])
#  N[i] ~ Poisson(lambda[i])
#Replicated counts are simulated as usual under a binomial observation model, and detection probability is allowed to vary by one site and one observational covariate: respectively quadrat forest cover, which is real data in the BerneseOberland data set, and wind-speed, which is invented data.
#Counts at each site (i) and for each occasion (j) are then produced according to the following model:
#  p[i,j] <- plogis(alpha0 + alpha1 * forest[i] + alpha2 * wind[i,j])
#  C[i,j] ~ Binomial(N[i], p[i,j])
#Finally, we assume that not each one of the 2500 quadrats is surveyed. Hence, we allow you to choose the number of quadrats that are surveyed and these will then be randomly placed into the landscape. We then assume that counts will only be available for these surveyed quadrats, i.e., counts from all non-surveyed quadrats will be NA'd out.

# truncN is a graphical parameter that truncates the z axis in some plots

BerneseOberland <- NULL # otherwise "no visible binding for global variable 'BerneseOberland'" when checked
data(BerneseOberland, envir = environment())

# Simulate spatial random field
set.seed(seeds[1])
s <- simExpCorrRF(variance = variance.RF, theta = theta.RF, show.plots=show.plots)

# SimulateNmix data with spatially correlated random effect in abundance
nsites<- 2500    # Number of sites (corresponding to the 50 by 50 grid)
# nsurveys<- nsurveys  # Number of replicate observations
y <- array(dim = c(nsites, nsurveys)) # Array for counts

# Ecological process
beta0 <- log(mean.lambda)
elev <- standardize(BerneseOberland$elevation)
forest<- standardize(BerneseOberland$forest)
loglam0 <- beta0 + beta[1] * elev + beta[2] * elev^2
loglam<- loglam0 + c(s$field)
lam0 <- exp(loglam0)
lam<- exp(loglam)

# Determine actual abundances as Poisson rv’s with parameter lam
N <- rpois(n = nsites, lambda = lam)
sum(N > 0) / nsites          # Finite-sample occupancy
Ntotal<- sum(N)

# Observation process
# Detection probability as linear function of forest and wind speed
alpha0 <- qlogis(mean.p)
wind<- matrix(rnorm(nsites*nsurveys), nrow = nsites, ncol = nsurveys)
p <- array(NA, dim = c(nsites, nsurveys))
for(j in 1:nsurveys){
  p[,j]<- plogis(alpha0 + alpha[1] * forest + alpha[2] * wind[,j])
}

# Go out and count things
for (j in 1:nsurveys){
  y[,j] <- rbinom(n = nsites, size = N, prob = p[,j])
}

# Select a sample of sites for surveys
set.seed(seeds[2])
surveyed.sites<- sort(sample(1:nsites, size = sample.size))

# Create the array of observed data by NA'ing out unsurveyed quadrats
yobs<- y    # Make a copy: the observed data
yobs[-surveyed.sites,] <- NA


# Minimal console output
true <- sum(N)
obs <- sum(apply(y, 1, max))
if(verbose) {
  cat("\n\nTrue total population size:", true)
  cat("\nTheoretically observed population size (sumMaxC) in 2500 quadrats:",obs)
  cat(paste("\nObserved population size in", sample.size,"surveyed quadrats:", sum(apply(yobs, 1, max), na.rm = TRUE)))
  cat("\nUnderestimation of abundance in total of 2500 quadrats:", round(100*(1-obs/true)), "%\n\n")
}
# Plot stuff
if(show.plots){
  # Restore graphical settings on exit ---------------------------
  oldpar <- par(mfrow = c(1,2), mar = c(5,8,5,2), cex.lab = 1.5, "cex.main")
  oldAsk <- devAskNewPage(ask = dev.interactive(orNone=TRUE))
  on.exit({par(oldpar); devAskNewPage(oldAsk)})
  # --------------------------------------------------------------

  tryPlot <- try( {
    # Plot lambda as function of covariates (excluding spatial field)
    plot(BerneseOberland$elevation, lam0, cex = 1, pch = 16,
        main = "Expected counts (lambda) without spatial field", xlab = "Elevation",
        ylab = "lambda excl. spatial field", frame = FALSE, col = rgb(0, 0, 0, 0.3))
    # Plot lambda as function of covariates (with spatial field)
    plot(BerneseOberland$elevation, lam, cex = 1, pch = 16,
        main = "Expected counts (lambda) with effect of spatial field", xlab = "Elevation",
        ylab = "lambda incl. spatial field", frame = FALSE, col = rgb(0, 0, 0, 0.3))

    # Plot detection as a function of the two covariates
    par(mfrow = c(1,2), cex.main = 1.5)
    plot(wind, p, ylim = c(0,1), cex = 1, main = "Detection (p) ~ Wind speed",
      frame = FALSE, col = rgb(0,0,0,0.3), pch = 16)
    noforest<- forest < -1.34
    points(wind[noforest,], p[noforest,], col = 'blue', pch = 16, cex = 1)
    legend('topright', 'blue: sites with no forest')
    plot(BerneseOberland$forest, apply(p, 1, mean), ylim = c(0,1), cex = 1,
        main = "Detection (p) ~ Forest cover", frame = FALSE, col = rgb(0,0,0,0.3), pch = 16)

    # Summary set of plots
    par(mfrow = c(2, 3), mar = c(2,2,4,6))
    r <- raster::rasterFromXYZ(data.frame(x = BerneseOberland$x, y = BerneseOberland$y,
        z = BerneseOberland$elevation))
    raster::plot(r, col = topo.colors(20), axes = FALSE, box = FALSE,
        main = "Elevation (metres)")
    r <- raster::rasterFromXYZ(data.frame(x = BerneseOberland$x, y = BerneseOberland$y,
        z = BerneseOberland$forest))
    raster::plot(r, col = topo.colors(20), axes = FALSE, box = FALSE, main = "Forest cover (%)")
    r <- raster::rasterFromXYZ(data.frame(x = s$gr[,1], y = s$gr[,2], z = c(s$field)))
    raster::plot(r, col = topo.colors(20), axes = FALSE, box = FALSE, main = "Spatial effect (neg.exp. corr.)")
    r <- raster::rasterFromXYZ(data.frame(x = BerneseOberland$x, y = BerneseOberland$y, z = N))
    raster::plot(r, col = topo.colors(20), axes = FALSE, box = FALSE,
        main = paste("Abundance (N, truncated at", truncN, ")"), zlim = c(0, truncN))
    r <- raster::rasterFromXYZ(data.frame(x = BerneseOberland$x, y = BerneseOberland$y,
        z = apply(p, 1, mean)))
    raster::plot(r, col = topo.colors(20), axes = FALSE, box = FALSE, 
        main = "Average detection probability")
    r <- raster::rasterFromXYZ(data.frame(x = BerneseOberland$x, y = BerneseOberland$y, z = apply(y, 1, max)))
    raster::plot(r, col = topo.colors(20), axes = FALSE, box = FALSE,
        main = paste("Max count (truncated at", truncN, ")\n with surveyed quadrats"), zlim = c(0, truncN))
    points(BerneseOberland$x[surveyed.sites], BerneseOberland$y[surveyed.sites],
        pch = 16, col = "red", cex = 0.8)
  }, silent = TRUE)
  if(inherits(tryPlot, "try-error"))
    tryPlotError(tryPlot)
}

# Output
return(list(
  # ------------------- arguments input ----------------------
  nsurveys = nsurveys, mean.lambda = mean.lambda, beta = beta, mean.p = mean.p,
  alpha0 = alpha0, alpha = alpha, sample.size = sample.size, variance.RF = variance.RF,
  theta.RF = theta.RF, seeds = seeds,
  # ------------------- from BerneseOberland ----------------------
  xcoord = BerneseOberland$x,
  ycoord = BerneseOberland$y,
  elevation = BerneseOberland$elevation,
  forest = BerneseOberland$forest,
  elevationS = elev, forestS = forest,
  # ------------------- generated variables ------------------------
  wind = wind,
  field = s$field,
  beta0 = beta0,
  lam = lam,
  N = N,
  Ntotal = Ntotal,
  p = p,
  y = y,
  surveyed.sites = surveyed.sites,
  yobs = yobs))
} # End of function definition
