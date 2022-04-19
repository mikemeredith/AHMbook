# Functions for the book Applied Hierarchical Modeling in Ecology (AHM)
# Marc Kery & Andy Royle, Academic Press, 2016.

# simOcc - AHM1 section 10.5 p577

# Function to simulate data for static occupancy models under wide range of conditions
#   (introduced in AHM1 Section 10.5)
simOcc <- function(M = 267, J = 3, mean.occupancy = 0.6,
    beta1 = -2, beta2 = 2, beta3 = 1, mean.detection = 0.3, time.effects = c(-1, 1),
    alpha1 = -1, alpha2 = -3, alpha3 = 0, sd.lp = 0.5, b = 2, show.plots = TRUE){
#
# Written by Marc Kery, 21 March 2015
#
# Function to simulate occupancy measurements replicated at M sites
# during J occasions.
# Population closure is assumed for each site.
# Expected occurrence may be affected by elevation (elev),
# forest cover (forest) and their interaction.
# Expected detection probability may be affected by elevation,
# wind speed (wind) and their interaction.
# Function arguments:
#     M: Number of spatial replicates (sites)
#     J: Number of temporal replicates (occasions)
#     mean.occupancy: Mean occurrence at value 0 of occurrence covariates
#     beta1: Main effect of elevation on occurrence
#     beta2: Main effect of forest cover on occurrence
#     beta3: Interaction effect on occurrence of elevation and forest cover
#     mean.detection: Mean detection prob. at value 0 of detection covariates
#     time.effects (on logit scale): bounds for uniform distribution from
#     which time effects gamma will be drawn
#     alpha1: Main effect of elevation on detection probability
#     alpha2: Main effect of wind speed on detection probability
#     alpha3: Interaction effect on detection of elevation and wind speed
#     sd.lp: standard deviation of random site effects (on logit scale)
#     b: constant value of 'behavioural response' leading to 'trap-happiness'
#     (if b > 0) or 'trap shyness' (if b < 0)
#     show.plots: if TRUE, plots of the data will be displayed;
#        IMPORTANT: has to be set to FALSE if you are running simulations.
if(FALSE) x <- NULL # Fudge to stop R CMD check complaining about curve

# Checks and fixes for input data -----------------------------
M <- round(M[1])
J <- round(J[1])
stopifnotProbability(mean.occupancy)
stopifnotProbability(mean.detection)
stopifNegative(sd.lp)
# --------------------------------------------

# Create some data structures: observed data and 3 versions of p matrix
y <- p <- p0 <- p1 <- array(NA, dim = c(M,J))  # Create data structures

# Create 2 site covariates (elev, forest) and 1 obs. covar. (wind)
elev <- runif(n = M, -1, 1)                         # Scaled elevation
forest <- runif(n = M, -1, 1)                       # Scaled forest cover
wind <- array(runif(n = M*J, -1, 1), dim = c(M, J)) # Scaled wind speed

# Model for occurrence (presence/absence): simulate system state z
beta0 <- qlogis(mean.occupancy)            # Mean occurrence on link scale
psi <- plogis(beta0 + beta1*elev + beta2*forest + beta3*elev*forest)
z <- rbinom(n = M, size = 1, prob = psi)   # Realised occurrence (true state)

# Plots for system state moved to line

# Model for observations: simulate observations y, given system state z
alpha0 <- qlogis(mean.detection)        # mean detection on link scale
gamma <- runif(J, time.effects[1], time.effects[2]) # (fixed) time effects
eps <- rnorm(M, 0, sd.lp)               # Individual (random) effects


# Generate two full detection probability arrays
# p0: for no preceding capture, p1: for preceding capture event
for(j in 1:J){
   p0[,j] <- plogis(alpha0 + gamma[j] + alpha1*elev + alpha2*wind[,j] + alpha3*elev*wind[,j] + eps)     # p when not captured at occasion j-1
   p1[,j] <- plogis(alpha0 + gamma[j] + alpha1*elev + alpha2*wind[,j] + alpha3*elev*wind[,j] + eps + b) # p when captured at occasion j-1
}

# Generate detection/nondetection data: the measurements of presence/absence
# For the first capture occasion (no behavioural response possible)
p[,1] <- p0[,1]      # Write the detection probability matrix p
y[,1] <- rbinom(n = M, size = z, prob = p0[,1])
# y[,1] <- rbinom(n = M, size = 1, prob = z * p0[,1])  # SAME

# Later capture occasions (potentially with contribution of b)
for (j in 2:J){
   for(i in 1:M){
      p[i,j] <- (1-y[i,(j-1)])*p0[i,j] + y[i,(j-1)] * p1[i,j] # which p?
      y[i,j] <- rbinom(n = 1, size = 1, prob = z[i] * p[i,j])
   }
}

# True and observed measures of 'distribution'
sumZ <- sum(z)                     # Total occurrence (all sites)
sumZ.obs <- sum(apply(y,1,max))    # Observed number of occ sites
psi.fs.true <- sum(z) / M          # True proportion of occ. sites in sample
psi.fs.obs <- mean(apply(y,1,max)) # Observed proportion of occ. sites in sample

if(show.plots){
  # Restore graphical settings on exit -------------------------
  oldpar <- par("mfrow", "cex.main", "cex.lab", "mar")
  oldAsk <- devAskNewPage(ask = dev.interactive(orNone=TRUE))
  on.exit({par(oldpar); devAskNewPage(oldAsk)})
  # ------------------------------------------------------------

  tryPlot <- try( {
    # Plots for system state
    par(mfrow = c(2, 2), cex.main = 1)
    curve(plogis(beta0 + beta1*x), -1, 1, col = "red", frame.plot = FALSE, ylim = c(0, 1), xlab = "Elevation", ylab = "psi", lwd = 2)
    plot(elev, psi, frame.plot = FALSE, ylim = c(0, 1), xlab = "Elevation", ylab = "")
    curve(plogis(beta0 + beta2*x), -1, 1, col = "red", frame.plot = FALSE, ylim = c(0, 1), xlab = "Forest cover", ylab = "psi", lwd = 2)
    plot(forest, psi, frame.plot = FALSE, ylim = c(0, 1), xlab = "Forest cover", ylab = "")

    # Plots for observation process
    par(mfrow = c(2, 3), cex.main = 1.2, cex.lab = 1.5, mar = c(5,5,3,2))
    # Plots for elevation, time, 'heterogeneity', and 'behavioural response'
    # Plots for elevation and time
    curve(plogis(alpha0 + alpha1*x), -1, 1, col = "red", frame.plot = FALSE, ylim = c(0, 1),
      xlab = "Elevation", ylab = "Expected detection (p)", lwd = 2,
      main = "Effects of elev and time")
    for(j in 1:J){
      curve(plogis(alpha0 + gamma[j] + alpha1*x),-1,1,lwd = 1, col="grey", add=TRUE)
    }
    # Plots for elevation and 'heterogeneity'
    curve(plogis(alpha0 + alpha1*x), -1, 1, col = "red", frame.plot = FALSE, ylim = c(0, 1),
      xlab = "Elevation", ylab = "Expected detection (p)", lwd = 2,
      main = "Effects of elev and site heterogeneity")
    for(i in 1:M){
      curve(plogis(alpha0 + eps[i] + alpha1*x),-1,1,lwd = 1, col="grey", add=T)
    }
    curve(plogis(alpha0 + alpha1*x), -1, 1, col = "red", lwd = 2, add = TRUE)

    # Plot for elevation and 'behavioural response'
    p0plot <- p0
    p1plot <- p1   ;   p1plot[,1] <- NA
    for(j in 2:J){
       p0plot[,j] <- p0plot[,j] / (1 - y[,(j-1)])   # NA out some
       p1plot[,j] <- p1plot[,j] / y[,(j-1)]       # NA out some
    }
    matplot(elev, p0plot, xlab = "Elevation", ylab = "Detection (p)",
      main = "p ~ elevation at actual wind speed \n(red/blue - following/not following det.)",
      pch = 1, ylim = c(0,1), col = "blue", frame.plot = FALSE)
    if(sum(is.finite(p1plot)) > 0)
      matplot(elev, p1plot, pch = 16, col = "red", add = TRUE)


    # Plots for wind speed, time, 'heterogeneity', and 'behavioural response'
    # Plots for elevation and time
    curve(plogis(alpha0 + alpha2*x), -1, 1, col = "red", frame.plot = FALSE, ylim = c(0, 1),
      xlab = "Wind speed", ylab = "Expected detection (p)", lwd = 2,
      main = "Effects of wind and time")
    for(j in 1:J){
      curve(plogis(alpha0 + gamma[j] + alpha2*x),-1,1,lwd = 1, col="grey", add=TRUE)
    }
    # Plots for wind speed and 'heterogeneity'
    curve(plogis(alpha0 + alpha2*x), -1, 1, col = "red", frame.plot = FALSE, ylim = c(0, 1),
      xlab = "Wind speed", ylab = "Expected detection (p)", lwd = 2,
      main = "Effects of wind and site heterogeneity")
    for(i in 1:M){
      curve(plogis(alpha0 + eps[i] + alpha2*x),-1,1,lwd = 1, col="grey", add=TRUE)
    }
    curve(plogis(alpha0 + alpha2*x), -1, 1, col = "red", lwd = 2, add = TRUE)

    # Plot for wind speed and 'behavioural response'
    matplot(wind, p0plot, xlab = "Wind speed", ylab = "Detection (p)",
      main = "p ~ elevation at actual elevation \n(red/blue - following/not following det.)",
      pch = 1, ylim = c(0,1), col = "blue", frame.plot = FALSE)
    if(sum(is.finite(p1plot)) > 0)
      matplot(wind, p1plot, pch = 16, col = "red", add = TRUE)
  }, silent = TRUE)
  if(inherits(tryPlot, "try-error"))
    tryPlotError(tryPlot)
}

# Output
return(list(M = M, J = J, mean.occupancy = mean.occupancy, beta0 = beta0,
    beta1 = beta1, beta2 = beta2, beta3 = beta3, mean.detection = mean.detection,
    time.effects = time.effects, gamma = gamma, alpha0 = alpha0, alpha1 = alpha1,
    alpha2 = alpha2, alpha3 = alpha3, sd.lp = sd.lp, eps = eps, b = b,
    elev = elev, forest = forest, wind = wind, psi = psi, z = z,
    p = p, p0 = p0, p1 = p1, y = y, sumZ = sumZ, sumZ.obs = sumZ.obs,
    psi.fs.true = psi.fs.true, psi.fs.obs = psi.fs.obs))
}

