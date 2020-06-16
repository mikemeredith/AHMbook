
# Another function for Chapter 15 in AHM2

simDemoDynocc<- function(nsites = 100, nyears = 10, nvisits = 5, psi1 = 0.6,
    range.phi = c(0.2, 0.9), range.r = c(0, 0.4), range.p = c(0.1, 0.9),
    show.plot=TRUE) {
  #
  # Function simulates data under a variant of the demographic occupancy
  #  (or 'local survival') model of Roth & Amrhein (J. Appl. Ecol., 2010).
  #  Data are simulated in an 'unconditional' manner, i.e., for each site from first to last year.
  #  All parameter can be made year-dependent by specification of a range,
  #  within which annual values will be drawn from uniform distributions.
  #
  # What the function arguments mean:
  #   psi1 = probability a territory is occupied at t=1
  #   nsites = number of territories
  #   nyears = number of study years
  #   nvisits = number of replicate visits per site and year
  #   range.phi = lower and upper limit of uniform distribution, from which
  #     annual local survival probability is drawn
  #   range.r = lower and upper limit of uniform distribution, from which
  #     annual recruitment probability is drawn
  #   range.p = lower and upper limit of uniform distribution, from which
  #     annual detection probability is drawn

  # Checks and fixes for input data -----------------------------
  nsites <- round(nsites[1])
  nyears <- round(nyears[1])
  nvisits <- round(nvisits[1])
  stopifnotProbability(psi1)
  stopifnotProbability(range.phi) # bounds
  stopifnotProbability(range.r) # bounds
  stopifnotProbability(range.p) # bounds
  # ----------------------------------------------------------------

  # Define true territory occupancy state matrix z
  z <- matrix(rep(NA, nyears*nsites), ncol=nyears)

  # Define the 3-dimensional matrix y that contains the observations
  y <- array(NA, dim = c(nsites, nvisits, nyears))

  # Simulate the annual local survival (nyears - 1 intervals)
  phi <- runif(nyears-1, min(range.phi), max(range.phi))

  # Simulate the annual colonization (nyears - 1 intervals)
  r <- runif(nyears-1, min(range.r), max(range.r))

  # Simulate the annual detection (includes year 1)
  p <- runif(nyears, min(range.p), max(range.p))

  # Simulate true state z from t=1:nyears
  persistence <- new.colonization <- z  # Provide intermediate structures
  for(i in 1:nsites) {
    # Initial year (t=1)
    z[i,1] <- rbinom(1, 1, psi1)
    for(t in 2:nyears) {
      persistence[i,t] <- z[i,t-1] * phi[t-1] +
        z[i,t-1] * (1-phi[t-1]) * r[t-1]  # survival or a 'rescue process'
      new.colonization[i,t] <- (1-z[i,t-1]) * r[t-1]
      z[i,t] <- rbinom(1, 1, persistence[i,t] + new.colonization[i,t])
    }
  }

  # Observations from t=1:nyears
  for(i in 1:nsites) {
    for(t in 1:nyears) {
      for(j in 1:nvisits) {
        y[i,j,t] <- rbinom(1, 1, z[i,t] * p[t])
      }
    }
  }

  # Create vector with 'occasion of marking' (for observed data)
  obsz <- apply(y, c(1,3), max)
  f <- suppressWarnings(apply(obsz, 1, function(x) min(which(x!=0))))
  f[f == 'Inf'] <- nyears

  # Derived quantities
  nocc.true <- apply(z, 2, sum)   # True ...
  nocc.obs <- apply(obsz, 2, sum) # ... and observed number of pairs

  if(show.plot) {
    # Visualization by two graphs
    oldpar <- par(mfrow = c(1, 2), mar = c(5, 5, 4, 2), cex.lab = 1.5)
      on.exit(par(oldpar))
    tryPlot <- try( {
      plot(1, 0, type = 'n', ylim = c(0,1), frame = FALSE,
          xlab = "Year", ylab = "Probability", xlim = c(1, nyears), las = 1,
          main = 'Local survival, recruitment and detection', xaxt='n')
      axis(1, 1:nyears)
      lines(1:(nyears-1), phi, type = 'o', pch=16, lwd = 2, col = 4, lty=2)
      lines(1:(nyears-1), r, type = 'o', pch=16, lwd = 2, col = 2, lty=3)
      lines(1:nyears, p, type = 'o', pch=16, lwd = 2, col = 1)
      legend('top', c("survival", "recruitment", "detection"),
          lty=c(2,3,1), lwd=2, col=c(4,2,1), #pch=16,
          inset=c(0, -0.05), bty='n', xpd=NA, horiz=TRUE)

      plot(1:nyears, nocc.true, type = 'n', frame = FALSE, xlab = "Year",
          ylab = "Population size", xlim = c(1, nyears), ylim = c(0, nsites), las = 1,
          main = 'True and observed population size', xaxt='n')
      axis(1, 1:nyears)
      lines(1:nyears, nocc.true, type = 'o', pch=16, lwd = 2, col = 2, lty=1)
      lines(1:nyears, nocc.obs, type = 'o', pch=16, lwd = 2, col = 4, lty=2)
      legend('top', c("true", "observed"),
          lty=c(1,2), lwd=2, col=c(2,4),
          inset=c(0, -0.05), bty='n', xpd=NA, horiz=TRUE)
    }, silent = TRUE)
    if(inherits(tryPlot, "try-error"))
      tryPlotError(tryPlot)
  }

  # Return stuff
  return(list(
    # ----------- arguments supplied -----------------------
    psi1 = psi1, nsites = nsites, nyears = nyears, nvisits = nvisits,
    range.phi = range.phi, range.r = range.r, range.p = range.p,
    # ----------- generated values ---------------------------
    phi = phi,
    r = r,
    p = p,
    z = z,
    y = y,
    f = f,
    nocc.true = nocc.true,
    nocc.obs = nocc.obs))
}

