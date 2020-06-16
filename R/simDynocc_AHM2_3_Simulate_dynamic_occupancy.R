# AHM2 chapter 4

# Revised 4 Dec 2018, 1 March 2019

# ------------------ Start function definition ---------------------
simDynocc<- function(nsites = 250, nyears = 10, nsurveys = 3, year.of.impact = NA,
    mean.psi1 = 0.4, beta.Xpsi1 = 0,
    range.phi = c(0.5, 1), sd.lphi.site = 0, impact.phi = 0, beta.Xphi = 0,
    range.gamma = c(0, 0.5), sd.lgamma.site = 0, impact.gamma = 0, beta.Xgamma = 0,
    sd.lphi.lgamma.site = 0,
    range.p = c(0.1, 0.9), beta.Xp = 0,
    range.beta1.survey = c(0, 0), range.beta2.survey = c(0, 0), trend.sd.site = c(0, 0),
    trend.sd.survey = c(0, 0), trend.sd.site.survey = c(0, 0), show.plots = TRUE) {
  #
  # Written by Marc Kery, 4 Dec 2014 - 4 December 2018
  #
  # Function to simulate detection/nondetection data under a general
  #     dynamic site-occ model, including:
  #   * annual variation in the probabilities of site persistence and colonization
  #     and detection is specified by the bounds of a uniform distribution.
  #   * one covariate is allowed to affect a parameter: a site covariate for psi1,
  #     a site-by-year covariate for phi and gamma and an
  #     observational covariate for p
  #   * Additional detection heterogeneity at the site, survey (= occasion),
  #     or the site-by-survey level, with the possibility of a temporal trend in
  #     these types of heterogeneity over the years.
  #        For instance, an annual trend in detection heterogeneity at the site
  #     or the survey level is specified by different first and second values,
  #     which correspond to the heterogeneity in the first and the last year,
  #     with a linear trend interpolated for the intervening years.
  #        As an example, trend.sd.site = c(0, 1) will result in
  #     a linear trend in the magnitude of site-level heterogeneity in detection
  #     from 0 in the first year to 1 in the last year.
  #   * Additional detection heterogeneity that varies over the season
  #     (= among occasions) according to a linear and quadratic occasion effect
  #     (e.g., to model the phenology of an insect species).
  #   * Data simulation under a BACI (before-after-control-impact) design is
  #     possible,where an event happens in a specified year and
  #     reduces phi or gamma by a stated percentage (only reductions possible)
  #
  # Function arguments:
  # -------------------
  # ** Sample size arguments **
  # nsites:  Number of sites
  # nyears:  Number of years (or 'seasons')
  # nsurveys:  Number of replicate surveys (= occasions) within a year
  #
  # ** Arguments to set intercepts of regressions **
  # mean.psi1 - average occupancy probability in first year
  # range.p - bounds of uniform distribution from which annual p drawn
  # range.psi and range.gamma - same for survival and colonization probability
  # -------------------
  #
  # ** Arguments to set slopes of regressions **
  # beta.Xpsi1, beta.Xphi, beta.Xgamma, beta.Xp - covariate coefficient in
  #      prob. of initial occupancy, persistence, colonization and detection.
  # -------------------
  #
  # *** Args.for detection heterogeneity among sites, occasions and site/survey
  # trend.sd.site: sd of normal distribution to model logit-normal noise in p
  #      at the site level in the first and the last year of the simulation,
  #      with value of intervening years interpolated.
  # trend.sd.survey: sd of normal distribution to model logit-normal noise in p
  #      ONLY at the rep = occasion = 'survey' level, in the first and
  #      the last year, with interpolation for intervening years
  # trend.sd.site.survey: sd of normal distribution to model logit-normal noise
  #      in p at the site/year/rep level, in the first and the
  #      last year, with interpolation in between
  # For these arguments, if the two values in the range are the
  #      same, a constant value is assumed over time, while if they are different,
  #      a linear trend is assumed over time.
  # range.beta1.occasion is the range of the annual variation in the linear effect
  #     of the survey occasion (e.g., of month 1-12 when nsurveys = 12)
  #     on detection (= product of availability and perceptibility)
  # range.beta2.occasion is the same for the quadratic effect of survey occasion
  # -------------------
  #
  # ** Arguments for the BACI design **
  # year.of.impact: year in which an impact happens,
  #      which can then affect phi and gamma; NA means no impact
  # impact.phi: negative effect in percent on annual phi, ignored if no impact;
  #      (e.g., impact.phi = 20 means a 20% reduction in phi)
  # impact.gamma: negative effect in percent on annual gamma, ignored if no impact.
  # Note that for the BACI design, nyears must be greater than 2 and
  #      year.of.impact must not be equal to the first or the last year

  # Checks and fixes for input data -----------------------------
  nsites <- round(nsites[1])
  nyears <- round(nyears[1])
  nsurveys <- round(nsurveys[1])
  year.of.impact <- round(year.of.impact[1])
  if(!is.na(year.of.impact))
    stopifnotGreaterthan(nyears, 2)
  stopifnotBetween(year.of.impact, 2, nyears-1, allowNA=TRUE)
  stopifnotProbability(range.phi) # bounds
  stopifnotBetween(impact.phi, 0, 100)
  stopifnotProbability(range.gamma) # bounds
  stopifnotBetween(impact.gamma, 0, 100)
  stopifnotProbability(range.p) # bounds
  stopifnotLength(trend.sd.site, 2) # trend
  stopifNegative(trend.sd.site)
  stopifnotLength(trend.sd.survey, 2) # trend
  stopifNegative(trend.sd.survey)
  stopifnotLength(trend.sd.site.survey, 2) # trend
  stopifNegative(trend.sd.site.survey)
  # ----------------------------------------------------------------

  # Set up arrays needed
  site <- 1:nsites                        # Sites
  year <- 1:nyears                        # Years
  visit <- 1:nsurveys                     # Surveys (= months, visits, occasions)
  psi <- muZ <- z <- array(dim = c(nsites, nyears),
    dimnames = list(paste('Site', site, sep = ''), paste('Year', year, sep = ''))) # Occupancy, occurrence
  phi <- gamma <- array(NA, dim = c(nsites, (nyears-1)),
    dimnames = list(paste('Site', site, sep = ''), paste('Year', year[-nyears], sep = ''))) # Survival, colonisation
  y <- p <- array(NA, dim = c(nsites, nsurveys, nyears),
    dimnames = list(paste('Site', site, sep = ''), paste('Visit', visit, sep = ''),
      paste('Year', year, sep = '')))# Det. hist and p

  # Create covariates
  # Site covariate for psi1
  Xpsi1 <- runif(nsites, -2, 2)

  # Yearly-site covariates for phi and gamma
  Xphi <- array(runif(nsites*nyears, -2, 2), dim = c(nsites,nyears))
  Xgamma <- array(runif(nsites*nyears, -2, 2), dim = c(nsites,nyears))

  # Observational covariate for p
  Xp <- array(runif(nsites*nsurveys*nyears,-2,2),dim=c(nsites,nsurveys,nyears))

  # Create impact covariate for the BACI effect on phi and gamma
  impact  <- rep(0, (nyears-1))
  if(!is.na(year.of.impact))
    impact[year.of.impact:(nyears-1)] <- 1

  # (1) Simulate all parameter values
  # (a) State process parameters
  psi[,1] <- plogis(qlogis(mean.psi1) + beta.Xpsi1 * Xpsi1)
  mean.phi <- runif(n = nyears-1, min = min(range.phi), max = max(range.phi))
  mean.gamma <- runif(n = nyears-1, min = min(range.gamma), max = max(range.gamma))
  # new 2019-03-01: heterogeneity in phi and gamma across sites
  eps.lphi.site <- rnorm(n = nsites, mean = 0, sd = sd.lphi.site)
  eps.lgamma.site <- rnorm(n = nsites, mean = 0, sd = sd.lgamma.site)
  eps.lphi.lgamma.site <- rnorm(n = nsites, mean = 0, sd = sd.lphi.lgamma.site)

  # BACI effect on phi and gamma: negative effect on persistence/colonisation
  BACI.effect.phi <- (mean.phi * (impact.phi/100) * impact)
  BACI.effect.gamma <- (mean.gamma * (impact.gamma/100) * impact)
    # These will be zero if no BACI impact

  # Assemble effects of year, impact and covariates on phi and gamma
  for(t in 1:(nyears-1)){
    phi[,t] <- plogis(qlogis(mean.phi[t] - BACI.effect.phi[t]) +
      eps.lphi.site + eps.lphi.lgamma.site + beta.Xphi * Xphi[,t])
    gamma[,t] <- plogis(qlogis(mean.gamma[t] - BACI.effect.gamma[t]) +
      eps.lgamma.site + eps.lphi.lgamma.site + beta.Xgamma * Xgamma[,t])
  }

  # (b) Observation process parameters
  mean.p <- runif(n = nyears, min = min(range.p), max = max(range.p))
  beta1 <- runif(n = nyears, min = min(range.beta1.survey), max = max(range.beta1.survey))
  beta2 <- runif(n = nyears, min = min(range.beta2.survey), max = max(range.beta2.survey))

  # Next two allow incorporation of trend over time
  sd.site <- seq(from = trend.sd.site[1], to = trend.sd.site[2], length.out = nyears)
  sd.survey <- seq(from = trend.sd.survey[1], to = trend.sd.survey[2], length.out = nyears)
  sd.site.survey <- seq(from = trend.sd.site.survey[1], to = trend.sd.site.survey[2],
    length.out = nyears)

  # Define and fill the array of site/survey random effects
  eps3 <- array(dim = c(nsites, nsurveys, nyears))
  for(t in 1:nyears){
    eps3[,,t] <- matrix(rnorm(n = nsites*nsurveys, sd = sd.site.survey[t]), ncol = nsurveys)
  }

  for(i in 1:nsites){     # Sites
    for(t in 1:nyears){   # Years
      eps1 <- rnorm(n = nsites, sd = sd.site[t])   # Zero-mean site random eff.
      eps2 <- rnorm(n = nsurveys, sd = sd.survey[t]) # Zero-mean survey random eff.
      # ZM site.survey ranef.
      for(j in 1:nsurveys){ # Months
        p[i,j,t] <- plogis(qlogis(mean.p[t]) + beta.Xp*Xp[i,j,t] +
              eps1[i] + eps2[j] + eps3[i,j,t] +
        beta1[t] * (j - (nsurveys/2)) + beta2[t] * (j - (nsurveys/2))^2)
      }
    }
  }

  # (2) Simulate the true system dynamics (state process)
  # First year
  z[,1] <- rbinom(nsites, 1, psi[,1])   # Initial occurrence state
  # Years 2:nyears
  for(i in 1:nsites){                   # Loop over sites
    for(t in 2:nyears){                 # Loop over years
      muZ[i,t] <- z[i, t-1]*phi[i,t-1] + (1-z[i, t-1])*gamma[i,t-1]
      z[i,t] <- rbinom(1, 1, muZ[i,t])
    }
  }

  # (3) Simulate observation process to get the observed data
  for(i in 1:nsites){                     # Loop over sites
    for(t in 1:nyears){                   # Loop over years
      for(j in 1:nsurveys){               # Loop over replicates
        prob <- z[i,t] * p[i,j,t]        # zero out p for unoccupied sites
        y[i,j,t] <- rbinom(1, 1, prob)
      }
     }
  }

  # (4) Compute annual population occupancy
  for(i in 1:nsites){
    for (t in 2:nyears){
      psi[i,t] <- psi[i,t-1]*phi[i,t-1] + (1-psi[i,t-1])*gamma[i,t-1]
    }
  }
  n.occ <- apply(z, 2, sum)                         # True number of occupied sites
  n.occ.ever <- sum(apply(z, 1, max))               # True number of occupied sites ever
  zobs <- apply(y, c(1,3), max)                     # Observed presence-absence matrix
  n.occ.obs <- apply(zobs, 2, sum)                  # Observed number of occupied sites
  n.occ.ever.obs <- sum(apply(zobs, 1, max))        # Obs. number of occupied sites ever
  psi.fs <- apply(z, 2, mean)                       # Finite-sample occupancy proportion
  mean.psi <- apply(psi, 2, mean)                   # Average psi over sites
  psi.app <- apply(apply(y, c(1,3), max), 2, mean)  # Apparent occupancy (finite sample)

  # Compute average product of availability and detection in each occasion
  # (ignoring the other terms in the model for detection)
  p.occasion <- array(NA, dim = c(nsurveys, nyears))
  for(t in 1:nyears){   # Years
    p.occasion[,t] <- plogis(qlogis(mean.p[t]) + beta1[t] * (visit - (nsurveys/2)) +
        beta2[t] * (visit - (nsurveys/2))^2)
  }

  # Compute annual average of phi, gamma and p
  avg.phi <- colMeans(phi)
  avg.gamma <- colMeans(gamma)
  avg.p <- colMeans(apply(p, c(1,3), mean))

  # (5) Plots of stuff
  if(show.plots){
    # Restore graphical settings on exit
    oldpar <- par("mfrow", "mar", "cex.lab", "cex.axis")
    oldAsk <- devAskNewPage(ask = dev.interactive(orNone=TRUE))
    on.exit({par(oldpar); devAskNewPage(oldAsk)})

    tryPlot <- try( {
      # ------ Plot A ---------
      par(mfrow = c(2, 2), mar = c(5,5,5,3), cex.lab = 1.2, cex.axis = 1.2)
      # Get predicted covariate relationships and plot them in single graph
      pred.cov <- seq(-2, 2, length.out = 100)
      psi.pred <- plogis(qlogis(mean.psi1) + beta.Xpsi1 * pred.cov)
      phi.pred <- plogis(mean(qlogis(mean.phi)) + beta.Xphi * pred.cov)
      gamma.pred <- plogis(mean(qlogis(mean.gamma)) + beta.Xgamma * pred.cov)
      p.pred <- plogis(mean(qlogis(mean.p)) + beta.Xp * pred.cov)

      plot(pred.cov, psi.pred, type = 'n', xlim = c(-2, 2), ylim = c(0,1),
        main = 'Covariate relationships', xlab = 'Covariate value',
        ylab = 'Predicted probability', frame = FALSE, las=1)
      lines(pred.cov, psi.pred, type = 'l', col = 3, lwd = 3, lty=1)
      lines(pred.cov, phi.pred, type = 'l', col = 4, lwd = 3, lty=2)
      lines(pred.cov, gamma.pred, type = 'l', col = 1, lwd = 3, lty=3)
      lines(pred.cov, p.pred, type = 'l', col = 2, lwd = 1, lty=1)
      legend('top', legend = c('psi1', 'phi', 'gamma', 'p'),
        col = c(3,4,1,2), lty = c(1,2,3,1), lwd = c(3,3,3,1),
        inset=c(0, -0.15), bty='n', xpd=NA, horiz=TRUE)

      # Within-season pattern of detection (= product of availability and detection)
      # (ignoring the other terms in the model for detection)
      matplot(visit, p.occasion, type = 'l', lty = 1, lwd = 2,
        main = 'Within-season pattern in p over the years \n(only occasion terms)',
        xlab = 'Survey', ylab = 'Detection probability',
        ylim = c(0,1), frame = FALSE, xaxt='n')
      axis(1, at=1:nsurveys)

      # Histogram of detection
      hist(p, col = 'lightgrey', xlim = c(0,1), main = 'Detection probability p')

      # Plot realised and apparent proportion of occupied sites
      plot(year, avg.p, type = "n", xlab = "Year", ylab = "Probability",
        xlim = c(1,nyears), ylim = c(0,1), frame.plot = FALSE, las = 1, xaxt='n',
        main = 'True occupancy (finite-sample), \nobserved occupancy (prop. sites occupied) and average p')
      axis(1, 1:nyears)
      lines(year, apply(z, 2, mean), type = "l", col = 2, lwd = 2, lty = 1)
      lines(year, psi.app, type = "l", col = 1, lwd = 2, lty=2)
      lines(year, avg.p , type = "l", col = 2, lwd = 2, lty = 3)
      if(!is.na(year.of.impact)) {
        abline(v=year.of.impact+0.5, col='grey', lwd=2)
        text(year.of.impact+0.5, 0, "impact", adj=c(0.5, 0.5))#pos=1, offset=0)
      }
      legend('top', legend = c('True psi', 'Observed psi', 'Detection'),
        col = c(2,1,2), lty = c(1,2,3), lwd = 2,
        inset=c(0, -0.15), bty='n', xpd=NA, horiz=TRUE)

      # ------ Plot B ---------
      # Plot of population sizes (ever occupied, occupied per year, true and observed)
      # And comparison with annual vals of colonisation, persistence and detection
      par(mfrow = c(1, 2), mar = c(5,5,5,3), cex.lab = 1.2, cex.axis = 1.2)

      # Annual average of colonisation, persistence and detection
      plot(1:(nyears-1), avg.gamma, type = "n", xlab = "Year or Yearly interval",
        ylab = "Probability", xlim = c(0.5, nyears), ylim = c(0,1),
        las = 1, xaxt='n', frame.plot = FALSE,
        main = 'Average annual persistence,\ncolonization, and detection')
      axis(1, at=1:nyears)
      lines(1:(nyears-1), avg.phi, type = "o", pch=16, col = 4, lwd = 2, lty=3)
      lines(1:(nyears-1), avg.gamma, type = "o", pch=16, col = 1, lwd = 2, lty=2)
      lines(1:nyears, avg.p, type = "o", pch=16, col = 2, lwd = 2)
      if(!is.na(year.of.impact)) {
        abline(v=year.of.impact-0.5, col='grey', lwd=2)
        text(year.of.impact-0.5, 0, "impact", pos=1, offset=0)
      }
      legend('top', c("phi", "gamma", "p"),
        lty=c(3,2,1), lwd=2, col=c(4,1,2),
        inset=c(0, -0.05), bty='n', xpd=NA, horiz=TRUE)

      # True and observed number of occupied sites per year and overall (ever)
      plot(1, 1, type = "n", xlab = "Year",
        ylab = "Number of sites", xlim = c(0.5, nyears+1.5), xaxt='n',
        ylim = range(c(0, n.occ.ever, n.occ.obs)), frame.plot = FALSE,
        las = 1, main = 'True and obs. number of occupied sites \n per year and for all years combined (ever)')
      axis(1, at=1:nyears)
      end <- nyears/2 + 0.5
      mid <- mean(c(0.5, end))
      segments(0.5, n.occ.ever, end, n.occ.ever, lwd=2, col=2)
      text(mid, n.occ.ever, "True ever", pos=3, xpd=TRUE)
      segments(0.5, n.occ.ever.obs, end, n.occ.ever.obs, lwd=2, lty=2)
      text(mid, n.occ.ever.obs, "Observed ever", pos=1, xpd=TRUE)
      points(1:nyears, n.occ, type = "b", col = 2, pch = 16, cex = 1.5, lwd=3)
      points(1:nyears, n.occ.obs, type = "b", col = 1, pch=20, lty=2, cex = 1.5, lwd=3)
      if(!is.na(year.of.impact)) {
        abline(v=year.of.impact+0.5, col='grey', lwd=2)
        text(year.of.impact+0.5, 0, "impact", pos=1, offset=0)
      }
      legend('topright', legend = c('True annual', 'Obs. annual'),
        col = c(2,1), pch = c(16, 16), lty = c(1,1), pt.cex=1.5,
        lwd = 3, bty = 'n')
    }, silent = TRUE)
    if(inherits(tryPlot, "try-error"))
      tryPlotError(tryPlot)
  }

  # Return data
  return(list(nsites=nsites, nyears=nyears, nsurveys=nsurveys,year.of.impact = year.of.impact, impact = impact, mean.psi1=mean.psi1, beta.Xpsi1=beta.Xpsi1,
  range.phi=range.phi,
  sd.lphi.site = sd.lphi.site, impact.phi = impact.phi, beta.Xphi=beta.Xphi, BACI.effect.phi = BACI.effect.phi,
  range.gamma=range.gamma, sd.lgamma.site = sd.lgamma.site, impact.gamma = impact.gamma, beta.Xgamma=beta.Xgamma, sd.lphi.lgamma.site = sd.lphi.lgamma.site,
  BACI.effect.gamma = BACI.effect.gamma, range.p=range.p, beta.Xp=beta.Xp, trend.sd.site=trend.sd.site, trend.sd.survey=trend.sd.survey, range.beta1.survey = range.beta1.survey, range.beta2.survey = range.beta2.survey, beta1 = beta1, beta2 = beta2, p.occasion = p.occasion,sd.site=sd.site, sd.survey=sd.survey, mean.phi=mean.phi, mean.gamma=mean.gamma, mean.p=mean.p, avg.phi=avg.phi, avg.gamma=avg.gamma, avg.p=avg.p, psi=psi, mean.psi=mean.psi, n.occ = n.occ, n.occ.ever = n.occ.ever, n.occ.obs = n.occ.obs, n.occ.ever.obs = n.occ.ever.obs,
  psi.fs = psi.fs, psi.app=psi.app, z=z, phi=phi, gamma=gamma, p=p, y = y, Xpsi1 = Xpsi1, Xphi = Xphi, Xgamma = Xgamma, Xp = Xp,
  eps.lphi.site = eps.lphi.site, eps.lgamma.site = eps.lgamma.site, eps.lphi.lgamma.site = eps.lphi.lgamma.site,
  eps1 = eps1, eps2 = eps2, eps3 = eps3))
} # ------------------ End function definition ---------------------
