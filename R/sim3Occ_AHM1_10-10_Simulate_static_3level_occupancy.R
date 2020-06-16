# Functions for the book Applied Hierarchical Modeling in Ecology (AHM)
# Marc Kery & Andy Royle, Academic Press, 2016.

# sim3Occ - AHM1 section 10.10 p604

# Function to simulate data for static 3-level occupancy models
#   (introduced in AHM1 Section 10.10)
sim3Occ <- function(nunits = 100, nsubunits = 5, nreps = 3, mean.psi = 0.8, beta.Xpsi = 1,
  sd.logit.psi = 0, mean.theta = 0.6, theta.time.range = c(-1, 1), beta.Xtheta = 1,
  sd.logit.theta = 0, mean.p = 0.4, p.time.range = c(-2,2), beta.Xp = -1, sd.logit.p = 0,
  show.plot = TRUE, verbose = TRUE) {
#
# Function generates 3-level occupancy data
#   with possibility of site-specific random variation at every level,
#   "time effects" at the middle and the lower levels and
#   effects of one distinct covariate at each level.
#
# Written by Marc Kery, 2014
#
# Function arguments:
# nunits: Number of main units (large quadrats)
# nsubunits: Number of subunits (nested subsamples within
#    each main unit
# nreps: Number of rep surveys in every subunit
# mean.psi: Mean large-scale, unit-level occupancy probability (psi)
# beta.Xpsi: effect on psi of covariate A (at main unit level)
# sd.logit.psi: SD of logit(psi), unstructured site variation in psi
# mean.theta: Mean small-scale (subunit) occupancy probability (theta)
# theta.time.range: range of theta 'intercepts' for subunits
# beta.Xtheta: effect on theta of covariate B (at subunit level)
# sd.logit.theta: SD of logit(theta), unstructured site variation in theta
# mean.p: Mean per-survey detection probability
# p.time.range: range of p 'intercepts' for replicates
# beta.Xp: effect on p of covariate C (unit by subunit by replicate)
# sd.logit.p: SD of logit(p)
if(FALSE) x <- NULL # Fudge to stop R CMD check complaining.

  # Checks and fixes for input data -----------------------------
  nunits <- round(nunits[1])
  nsubunits <- round(nsubunits[1])
  nreps <- round(nreps[1])
  stopifnotProbability(mean.psi)
  stopifNegative(sd.logit.psi)
  stopifnotProbability(mean.theta)
  stopifNegative(sd.logit.theta)
  stopifnotProbability(mean.p)
  stopifNegative(sd.logit.p)
  # ----------------------------------------------------------------

# Create data structures
z <- psi <- array(NA, dim = nunits)  # Unit occurrence
a <- theta <- array(NA, dim = c(nunits, nsubunits)) # Subunit
y <- p <- array(NA, dim=c(nunits, nsubunits, nreps) ) # Rep

# Create standardised covariate values
covA <- as.numeric(array(runif(nunits, -2, 2), dim = nunits))
covB <- array(runif(nunits*nsubunits, -2, 2),
   dim = c(nunits, nsubunits))
covC <- array(runif(nunits*nsubunits*nreps, -2, 2),
   dim=c(nunits, nsubunits, nreps) )

# Simulate psi, theta and p and plot all
psi <- plogis(qlogis(mean.psi) + beta.Xpsi * covA + rnorm(nunits, 0, sd.logit.psi))
theta.time.effect <- runif(nsubunits, theta.time.range[1], theta.time.range[2])
p.time.effect <- runif(nreps, p.time.range[1], p.time.range[2])

for(j in 1:nsubunits){
   theta[,j] <- plogis(qlogis(mean.theta) + theta.time.effect[j]+ (beta.Xtheta*covB)[,j] + array(rnorm(nunits*nsubunits, 0, sd.logit.theta), dim = c(nunits, nsubunits))[,j])
   for(k in 1:nreps){
      p[,j,k] <- plogis(qlogis(mean.p) + p.time.effect[k] + (beta.Xp*covC)[,j,k]+ array(rnorm(nunits*nsubunits*nreps, 0,sd.logit.p),dim =c(nunits, nsubunits, nreps))[,j,k])
   }
}

# Visualisation of covariate relationships of psi, theta and p
if(show.plot) {
  op <- par(mfrow = c(1,3), mar = c(5,5,5,2), cex.lab = 1.5, cex.axis = 1.5) ; on.exit(par(op))
  tryPlot <- try( {
    plot(covA, psi, xlab = "Unit covariate A", ylab = "psi", ylim = c(0,1), main = "Large-scale occupancy probability (psi)", frame = FALSE)
    curve(plogis(qlogis(mean.psi) + beta.Xpsi * x), -2, 2, col = "red", lwd = 3, add = TRUE)
    plot(covB, theta, xlab = "Unit-subunit covariate B", ylab = "theta", ylim = c(0,1), main = "Small-scale occupancy probability/availability \n(theta) (red - time variation)", frame = FALSE)
    for(j in 1:nsubunits){
       curve(plogis(qlogis(mean.theta) + theta.time.effect[j] +
       beta.Xtheta * x), -2, 2, lwd = 2, col = "red", add = T)
    }
    plot(covC, p, xlab = "Unit-subunit-rep covariate C", ylab = "p", ylim = c(0,1), main = "Detection probability (p) \n (red - replicate variation)", frame = FALSE)
    for(k in 1:nreps){
       curve(plogis(qlogis(mean.p) + p.time.effect[k] +
       beta.Xp * x), -2, 2, lwd = 2, col = "red", add = T)
    }
  }, silent = TRUE)
  if(inherits(tryPlot, "try-error"))
    tryPlotError(tryPlot)
}
# Sample three nested Bernoulli distributions
# with probabilities psi, z*theta and a * p
for (i in 1:nunits) {
  z[i] <- rbinom(n = 1, size = 1, prob = psi[i])
  for (j in 1:nsubunits) {
    a[i, j] <- rbinom(n = 1, size = 1, prob = z[i] * theta[i,j])
    for (k in 1:nreps) {
      y[i,j,k] <- rbinom(n=1, size = 1, prob = a[i,j]*p[i,j,k])
    } # survey
  } # subunit
} # unit

sum.z <- sum(z)
sum.z.a <- sum(apply(a, 1, sum)>0)
obs.sum.z <- sum(apply(apply(y, c(1,2), max), 1, max))
if(verbose) {
  cat(" Occupied units:                           ", sum.z, "\n",
      "Units with >=1 occupied, surveyed subunit:", sum.z.a, "\n",
      "Observed number of occupied units:        ", obs.sum.z, "\n",
      "\n")
}
# Output
return(list(nunits = nunits, nsubunits = nsubunits, nreps = nreps, mean.psi = mean.psi, beta.Xpsi = beta.Xpsi, sd.logit.psi = sd.logit.psi, psi = psi, mean.theta = mean.theta, theta.time.range = theta.time.range, theta.time.effect = theta.time.effect, beta.Xtheta = beta.Xtheta, sd.logit.theta = sd.logit.theta, theta = theta, mean.p = mean.p, p.time.range = p.time.range, p.time.effect = p.time.effect, beta.Xp = beta.Xp, sd.logit.p = sd.logit.p, p = p, z = z, a = a, y = y, sum.z = sum.z, obs.sum.z = obs.sum.z, sum.z.a = sum.z.a, covA = covA, covB = covB, covC = covC))
}

