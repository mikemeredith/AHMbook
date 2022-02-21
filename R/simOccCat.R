# Functions for the book Applied Hierarchical Modeling in Ecology (AHM)
# Marc Kery & Andy Royle, Academic Press, 2016.

# SimOccCat = adaptation of simOcc (AHM1 section 10.5 p577) to allow for categorical covariates

# Function to simulate data for static occupancy models under wide range of conditions

simOccCat <- function(M = 267, J = 3, mean.occupancy = 0.6,
    beta1 = 0, beta2 = 0, beta3 = 0, mean.detection = 0.3, time.effects = c(0, 0),
    alpha1 = 0, alpha2 = 0, alpha3 = 0, sd.lp = 0, b = 0,
    nHab = 5, range.HAB = 2, nObs = 10, range.OBS = 4,   # new arguments
    show.plots = TRUE){
  #
  # The new arguments are:
  #   nHab: the number of categories for the site covariate
  #   range.HAB: controls the size of the effects for the categories of the site covariate
  #   nObs: the number of categories for the detection covariate
  #   range.OBS: controls the size of the effects for the categories of the detection covariate

  if(FALSE) x <- NULL # Fudge to stop R CMD check complaining about 'curve'

  # Checks and fixes for input data -----------------------------
  M <- round(M[1])
  J <- round(J[1])
  stopifnotProbability(mean.occupancy)
  stopifnotProbability(mean.detection)
  stopifNegative(sd.lp)
  # MORE TODO
  # --------------------------------------------

  # Create 2 continuous site covariates (elev, forest) and 1 obs. covar. (wind)
  elev <- runif(n = M, -1, 1)                         # Scaled elevation
  forest <- runif(n = M, -1, 1)                       # Scaled forest cover
  wind <- array(runif(n = M*J, -1, 1), dim = c(M, J)) # Scaled wind speed

  # Create categorical covariates with approximately equal sized categories
  HAB <- sample(nHab, M, replace = TRUE)
  OBSvec <- sample(nObs, M*J, replace = TRUE) # No constraint that observers only visit sites once

  # Create coefficients for HAB factor that sum to 0, calculate HAB effect
  coefHAB <- runif(nHab, 0, range.HAB)
  coefHAB <- coefHAB - mean(coefHAB)  # Now sums to zero
  HABeffect <- coefHAB[HAB]

  # Create coefficients for OBS factor that sum to 0, calculate OBS effect
  coefOBS <- runif(nObs, 0, range.OBS)
  coefOBS <- coefOBS - mean(coefOBS)  # Now sums to zero
  OBSeffect <- matrix(coefOBS[OBSvec], nrow=M, ncol=J)

  # Model for occurrence (presence/absence): simulate system state z
  beta0 <- qlogis(mean.occupancy)            # Mean occurrence on link scale
  psi <- plogis(beta0 + beta1*elev + beta2*forest + beta3*elev*forest + HABeffect)
  z <- rbinom(n = M, size = 1, prob = psi)   # Realised occurrence (true state)

  # Model for observations: simulate observations y, given system state z
  alpha0 <- qlogis(mean.detection)        # mean detection on link scale
  gamma <- runif(J, min(time.effects), max(time.effects)) # (fixed) time effects
  eps <- rnorm(M, 0, sd.lp)               # Site (random) effects


  # Generate detection probability array without behavioural effect
  # for(j in 1:J){
     # logit.p0[,j] <- alpha0 + gamma[j] + alpha1*elev + alpha2*wind[,j] + alpha3*elev*wind[,j] + eps + OBSeffect[,j]
  # }
  tmp <- alpha0 + alpha1*elev + alpha2*wind + alpha3*elev*wind + eps + OBSeffect
  logit.p0 <- sweep(tmp, 2, gamma, "+")

  # Generate detection/nondetection data: the measurements of presence/absence
  y <- p <- matrix(NA, M, J)
  # For the first capture occasion (no behavioural response possible)
  p[,1] <- plogis(logit.p0[,1])                    # 'p' is needed for the output
  y[,1] <- rbinom(n = M, size = z, prob = p[,1])
  # y[,1] <- rbinom(n = M, size = 1, prob = z * p0[,1])  # SAME

  # Later capture occasions (potentially with contribution of b)
  for (j in 2:J){
    p[, j] <- plogis(logit.p0[,j] + b*y[, j-1])
    y[, j] <- rbinom(n = M, size = z, prob = p[, j])
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
      p0plot <- plogis(logit.p0)
      p1plot <- plogis(logit.p0 + b)   ;   p1plot[,1] <- NA
         p0plot[,2:J] <- p0plot[,2:J] / (1 - y[,1:(j-1)])   # NA out some
         p1plot[,2:J] <- p1plot[,2:J] / y[,1:(j-1)]       # NA out some
      matplot(elev, p0plot, xlab = "Elevation", ylab = "Detection (p)",
        main = "p ~ elevation at actual wind speed \n(red/blue - following/not following det.)",
        pch = 1, ylim = c(0,1), col = "blue", frame.plot = FALSE)
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
      matplot(wind, p1plot, pch = 16, col = "red", add = TRUE)
    }, silent = TRUE)
    if(inherits(tryPlot, "try-error"))
      tryPlotError(tryPlot)
  }

  # Output
  return(list(
      # arguments input
      M = M, J = J, mean.occupancy = mean.occupancy, beta0 = beta0,
      beta1 = beta1, beta2 = beta2, beta3 = beta3, mean.detection = mean.detection,
      time.effects = time.effects, alpha0 = alpha0, alpha1 = alpha1,
      alpha2 = alpha2, alpha3 = alpha3, sd.lp = sd.lp, b = b,
      nHab = nHab, range.HAB = range.HAB, nObs = nObs, range.OBS = range.OBS,
      # Generated values
      gamma = gamma, eps = eps, elev = elev, forest = forest, wind = wind,
      HAB = HAB, OBS = matrix(OBSvec, M, J), coefHAB = coefHAB, coefOBS = coefOBS,
      psi = psi, z = z, p = p, p0 = plogis(logit.p0), p1 = plogis(logit.p0 + b), y = y,
      sumZ = sumZ, sumZ.obs = sumZ.obs, psi.fs.true = psi.fs.true, psi.fs.obs = psi.fs.obs))
}

