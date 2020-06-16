# Functions for the book Applied Hierarchical Modeling in Ecology (AHM)
# Marc Kery & Andy Royle, Academic Press, 2016.

# simOccttd - AHM1 section 10.12.1 p616

# Function simulates time-to-detection occupancy design data under model
# of Garrard et al. (Austral Ecology, 2008), also see Bornand et al. (MEE, 2014)
#   (introduced in AHM1 Section 10.12.1)

simOccttd <- function(M = 250, mean.psi = 0.4, mean.lambda = 0.3,
  beta1 = 1, alpha1 = -1, Tmax = 10, show.plot = TRUE, verbose = TRUE){
#
# Function simulates time-to-detection occupancy design data under model
# of Garrard et al. (Austral Ecology, 2008), also see Bornand et al. (MEE, 2014)
#
# Written by Marc Kery, 24 April 2015
#
# Function arguments:
# M: Number of sites
# mean.psi: intercept of occupancy probability
# mean.lambda: intercept of Poisson rate parameter
# beta1: slope of continuous covariate B on logit(psi)
# alpha1: slope of continuous covariate A on log(lambda)
# Tmax: maximum search time (in arbitrary units, which are same as response)
#   response will be censored at Tmax

# Checks and fixes for input data -----------------------------
M <- round(M[1])
stopifnotProbability(mean.psi)
stopifNegative(mean.lambda, allowZero=FALSE)
stopifNegative(Tmax, allowZero=FALSE)
# --------------------------------------------

# Generate covariate values
covA <- rnorm(M)
covB <- rnorm(M)

# Ecological process: Simulate occurrence z at each site
psi <- plogis(qlogis(mean.psi) + beta1 * covB)
(z <- rbinom(M, 1, psi))    # Realized occurrence at each site
if(verbose)
  cat("   Number of occupied sites ( among", M ,"):", sum(z), "\n")

# Observation process: Simulate time-to-detection (ttd) at each site
# Start without censoring
lambda <- exp(log(mean.lambda) + alpha1 * covA)
(ttd.temp <- rexp(M, lambda))   # ttd conditional on site occupied
# Now add two sources of censoringt censoring
ttd <- ttd.temp
ttd[z == 0] <- NA           # Censored if unoccupied
ttd[ttd.temp >= Tmax] <- NA # Censored if ttd >= Tmax

if(show.plot) {
  tryPlot <- try( {
    hist(ttd, breaks = length(ttd)/3, col = "gold",
        main = "Observed distribution of time to detection\n(censored cases (red line) excluded)",
        xlim = c(0, Tmax),
        xlab = "Measured time to detection")
    abline(v = Tmax, col = "red", lwd = 3)
  }, silent = TRUE)
  if(inherits(tryPlot, "try-error"))
    tryPlotError(tryPlot)
}

# Number of sites where detected
(n.obs <- sum(ttd < Tmax, na.rm = TRUE))
if(verbose)
  cat("   Number of sites at which detected:", n.obs, "\n")

# Calculate censoring indicator
d <- as.numeric(is.na(ttd))
if(verbose)
  cat("   Number of times censored:", sum(d), "\n")

# Output
return(list(M = M, mean.psi = mean.psi, mean.lambda = mean.lambda, beta1 = beta1, alpha1 = alpha1, Tmax = Tmax, covA = covA, covB = covB, psi = psi, lambda = lambda, z = z, ttd.temp = ttd.temp, ttd = ttd, d = d, sum.z = sum(z), n.obs = n.obs))
}



