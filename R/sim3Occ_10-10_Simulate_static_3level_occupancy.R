# Functions for the book Applied Hierarchical Modeling in Ecology (AHM)
# Marc Kéry & Andy Royle, Academic Press, 2016.

# sim3Occ - section 10.10 p604

# Function to simulate data for static 3-level occupancy models
#   (introduced in Section 10.10)
sim3Occ <- function(nunit = 100, nsubunit = 5, nrep = 3, mean.psi = 0.8, beta.Xpsi = 1, sd.logit.psi = 0, mean.theta = 0.6, theta.time.range = c(-1, 1), beta.Xtheta = 1, sd.logit.theta = 0, mean.p = 0.4, p.time.range = c(-2,2), beta.Xp = -1, sd.logit.p = 0){
#
# Function generates 3-level occupancy data
#   with possibility of site-specific random variation at every level,
#   "time effects" at the middle and the lower levels and
#   effects of one distinct covariate at each level.
#
# Written by Marc Kéry, 2014
#
# Function arguments:
# nunit: Number of main units (large quadrats)
# nsubunit: Number of subunits (nested subsamples within
#    each main unit
# nrep: Number of rep surveys in every subunit
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

# Create data structures
z <- psi <- array(NA, dim = nunit)  # Unit occurrence
a <- theta <- array(NA, dim = c(nunit, nsubunit)) # Subunit
y <- p <- array(NA, dim=c(nunit, nsubunit, nrep) ) # Rep

# Create standardised covariate values
covA <- as.numeric(array(runif(nunit, -2, 2), dim = nunit))
covB <- array(runif(nunit*nsubunit, -2, 2),
   dim = c(nunit, nsubunit))
covC <- array(runif(nunit*nsubunit*nrep, -2, 2),
   dim=c(nunit, nsubunit, nrep) )

# Simulate psi, theta and p and plot all
psi <- plogis(logit(mean.psi) + beta.Xpsi * covA + rnorm(nunit, 0, sd.logit.psi))
theta.time.effect <- runif(nsubunit, theta.time.range[1], theta.time.range[2])
p.time.effect <- runif(nrep, p.time.range[1], p.time.range[2])

for(j in 1:nsubunit){
   theta[,j] <- plogis(logit(mean.theta) + theta.time.effect[j]+ (beta.Xtheta*covB)[,j] + array(rnorm(nunit*nsubunit, 0, sd.logit.theta), dim = c(nunit, nsubunit))[,j])
   for(k in 1:nrep){
      p[,j,k] <- plogis(logit(mean.p) + p.time.effect[k] + (beta.Xp*covC)[,j,k]+ array(rnorm(nunit*nsubunit*nrep, 0,sd.logit.p),dim =c(nunit, nsubunit, nrep))[,j,k])
   }
}

# Visualisation of covariate relationships of psi, theta and p
par(mfrow = c(1,3), mar = c(5,5,5,2), cex.lab = 1.5, cex.axis = 1.5)
plot(covA, psi, xlab = "Unit covariate A", ylab = "psi", ylim = c(0,1), main = "Large-scale occupancy probability (psi)", frame = F)
curve(plogis(logit(mean.psi) + beta.Xpsi * x), -2, 2, col = "red", lwd = 3, add = TRUE)
plot(covB, theta, xlab = "Unit-subunit covariate B", ylab = "theta", ylim = c(0,1), main = "Small-scale occupancy probability/availability \n(theta) (red – time variation)", frame = F)
for(j in 1:nsubunit){
   curve(plogis(logit(mean.theta) + theta.time.effect[j] +
   beta.Xtheta * x), -2, 2, lwd = 2, col = "red", add = T)
}
plot(covC, p, xlab = "Unit-subunit-rep covariate C", ylab = "p", ylim = c(0,1), main = "Detection probability (p) \n (red – replicate variation)", frame = F)
for(k in 1:nrep){
   curve(plogis(logit(mean.p) + p.time.effect[k] +
   beta.Xp * x), -2, 2, lwd = 2, col = "red", add = T)
}

# Sample three nested Bernoulli distributions
# with probabilities psi, z*theta and a * p
for (i in 1:nunit) {
  z[i] <- rbinom(n = 1, size = 1, prob = psi[i])
  for (j in 1:nsubunit) {
    a[i, j] <- rbinom(n = 1, size = 1, prob = z[i] * theta[i,j])
    for (k in 1:nrep) {
      y[i,j,k] <- rbinom(n=1, size = 1, prob = a[i,j]*p[i,j,k])
    } # survey
  } # subunit
} # unit

sum.z <- sum(z)
sum.z.a <- sum(apply(a, 1, sum)>0)
obs.sum.z <- sum(apply(apply(y, c(1,2), max), 1, max))
cat(" Occupied units:                           ", sum.z, "\n",
    "Units with >=1 occupied, surveyed subunit:", sum.z.a, "\n",
    "Observed number of occupied units:        ", obs.sum.z, "\n",
    "\n")

# Output
return(list(nunit = nunit, nsubunit = nsubunit, nrep = nrep, mean.psi = mean.psi, beta.Xpsi = beta.Xpsi, sd.logit.psi = sd.logit.psi, psi = psi, mean.theta = mean.theta, theta.time.range = theta.time.range, theta.time.effect = theta.time.effect, beta.Xtheta = beta.Xtheta, sd.logit.theta = sd.logit.theta, theta = theta, mean.p = mean.p, p.time.range = p.time.range, p.time.effect = p.time.effect, beta.Xp = beta.Xp, sd.logit.p = sd.logit.p, p = p, z = z, a = a, y = y, sum.z = sum.z, obs.sum.z = obs.sum.z, sum.z.a = sum.z.a, covA = covA, covB = covB, covC = covC))
}

