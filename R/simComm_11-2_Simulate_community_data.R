# Functions for the book Applied Hierarchical Modeling in Ecology (AHM)
# Marc Kery & Andy Royle, Academic Press, 2016.

# simComm - section 11.2 p634

# Function to simulate community occupancy or community abundance data
#   with random species effects for psi/lambda and p (both including
#   effects of one covariate, 'habitat' for psi/lambda and 'wind speed' for p)
#   (introduced in Section 11.2)

simComm <- function(type=c("det/nondet", "counts"), nsite=30, nrep=3, nspec=100,
  mean.psi=0.25, sig.lpsi=1, mu.beta.lpsi=0, sig.beta.lpsi=0,
  mean.lambda=2, sig.loglam=1, mu.beta.loglam=1, sig.beta.loglam=1,
  mean.p=0.25, sig.lp=1, mu.beta.lp=0, sig.beta.lp=0, show.plot = TRUE) {
#
# Function simulates data from repeated sampling of a metacommunity
#   (or spatially structured community) according the model of
#   Dorazio & Royle (JASA, 2005) for type = "det/nondet" (this is the default)
#   or under the model of Yamaura et al. (2012) for type = "counts".
#
# Occupancy probability (psi) or expected abundance (lambda)
#   can be made dependent on a continuous site covariate 'habitat',
#   while detection probability can be made dependent an
#   observational covariate 'wind'.
# Both intercept and slope of the two log- or logistic regressions
#   (for occupancy or expected abundance, respectively, and for detection)
#   are simulated as draws from a normal distribution with
#   mean and standard deviation that can be selected using function arguments.
#
# Specifically, the data are simulated under the following linear models:
#
#  (1) for a type = "det/nondet" (i.e., community occupancy)
#  *********************************************************
#   (occupancy (psi) and detection (p) for site i, replicate j and species k)
#   psi[i,k] <- plogis(beta0[k] + beta1[k] * habitat[i]      # Occupancy
#   p[i,j,k] <- plogis(alpha0[k] + alpha1[k] * wind[i,j]     # Detection
#
#  (2) for a type = "counts" (i.e., community count)
#  ************************************************
#   (exp. abundance (lambda) and detection (p) for site i, rep. j and species k)
#   lambda[i,k] <- exp(beta0[k] + beta1[k] * habitat[i]      # E(N)
#   p[i,j,k] <- plogis(alpha0[k] + alpha1[k] * wind[i,j]     # Detection
#
#   Species-specific heterogeneity in intercepts and slopes is modelled
#   by up to four independent normal distributions (note: no correlation
#   between the intercepts as in Dorazio et al. (2006) or Kery & Royle (2008))
#
#  (1) for a type = "det/nondet" (i.e., community occupancy)
#  *********************************************************
#   beta0 ~ dnorm(qlogis(mean.psi), sig.lpsi)     # Mean and SD of normal distr.
#   beta1 ~ dnorm(mu.beta.lpsi, sig.beta.lpsi)
#   alpha0 ~ dnorm(qlogis(mean.p), sig.lp)
#   alpha1 ~ dnorm(mu.beta.lp, sig.beta.lp)
#
#  (2) for a type = "counts" (i.e., community count)
#  ************************************************
#   beta0 ~ dnorm(log(mean.lambda), sig.loglam)   # Mean and SD of normal distr.
#   beta1 ~ dnorm(mu.beta.loglam, sig.beta.loglam)
#   alpha0 ~ dnorm(qlogis(mean.p), sig.lp)
#   alpha1 ~ dnorm(mu.beta.lp, sig.beta.lp)


# Community occupancy model code partly based on code by Richard Chandler.
#
# Function arguments:
# *******************
# type: "det/nondet" or "counts"; hoose whether you want to
#    simulate detection/nondetection data or count data
# nsite: number of sites
# nrep: number of replicate samples (occasions or repeated measurements)
# nspec: total number of species in the area that is sampled by these sites
#    (regional species pool)
#
# mean.psi: community mean of occupancy probability over all species
#    in community (probability scale)
# sig.lpsi: community standard deviation of qlogis(occupancy probability intercept)
# mu.beta.lpsi: community mean of the effects of 'habitat'
#    covariate on occupancy probability
# sig.beta.lpsi: community standard deviation of the effects of
#    'habitat' covariate on occupancy probability
#
# mean.lambda: community mean of expected abundance over all species
#    in superpopulation
# sig.loglam: community standard deviation of log(lambda intercept)
# mu.beta.loglam: community mean of the effects of 'habitat' covariate
#    on log(lambda)
# sig.beta.loglam: community standard deviation of the effects
#    of 'habitat' covariate on expected abundance
#
# mean.p: community mean of detection probability over all species
#    in superpopulation (probability scale)
# sig.lp: community standard deviation of qlogis(detection probability intercept)
# mu.beta.lp: community mean of the effects of 'wind' covariate
#    on detection probability
# sig.beta.lp: community standard deviation of the effects of 'wind'covariate
#    on detection probability
# show.plot: choose whether to show plots or not. Set to FALSE when
#    using function in simulations.

# Code for simulating binary detection/nondetection data
#   (according to a community occupancy model)

if(FALSE) x <- NULL # A kludge to cope with 'curve's odd way of using 'x'

type <- match.arg(type)

if(show.plot){
  # Restore graphical settings on exit -------------------------
  oldpar <- par(no.readonly = TRUE)
  oldAsk <- devAskNewPage(ask = dev.interactive(orNone=TRUE))
  on.exit({par(oldpar); devAskNewPage(oldAsk)})
  # ------------------------------------------------------------
}

if(type=="det/nondet"){
# Prepare structures to hold data
  y.all <- y.obs <- p <- array(NA, c(nsite, nrep, nspec))
  dimnames(y.all) <- dimnames(y.obs) <- dimnames(p) <-
    list(paste("site", 1:nsite, sep=""),
    paste("rep", 1:nrep, sep=""),
    paste("sp", 1:nspec, sep=""))
  z <- psi <- matrix(NA, nsite, nspec)
  dimnames(z) <- dimnames(psi) <- list(paste("site", 1:nsite, sep=""),
    paste("sp", 1:nspec, sep=""))
  detected.at.all <- rep(NA, nspec)

  # Create covariates 'habitat' and 'wind'
  habitat <- sort(rnorm(nsite))     # Note 'habitat gradient' due to sorting
  wind <- matrix(rnorm(nsite * nrep), ncol=nrep)

  # Draw species-specific intercepts and slopes from their normal distributions
  # Build up linear predictors for occupancy and detection
  mu.lpsi <- ifelse(mean.psi == '1', 500, qlogis(mean.psi))
  mu.lp <- ifelse(mean.p == '1', 500, qlogis(mean.p))

  beta0 <- rnorm(nspec, mu.lpsi, sig.lpsi)           # occupancy intercept
  beta1 <- rnorm(nspec, mu.beta.lpsi, sig.beta.lpsi) # occupancy slope on habitat
  alpha0 <- rnorm(nspec, mu.lp, sig.lp)              # detection intercept
  alpha1 <- rnorm(nspec, mu.beta.lp, sig.beta.lp)    # detection slope on wind
  for(k in 1:nspec){
    psi[,k] <- plogis(beta0[k] + beta1[k] * habitat)
    for(j in 1:nrep){
      p[,j,k] <- plogis(alpha0[k] + alpha1[k] * wind[,j])
    }
  }

  # Distribute species over sites (simulate true state)
  for(k in 1:nspec){
    z[,k] <- rbinom(nsite, 1, psi[,k])
  }
  occurring.in.sample <- apply(z, 2, max) # Presence/absence at study sites

  # Measurement of presence/absence (simulate observation)
  for(k in 1:nspec) {
    for(i in 1:nsite){
      for(j in 1:nrep) {
        y.all[i,j,k] <- rbinom(1, z[i,k], p[i,j,k])
      }
    }
    # detected.at.all[k] <- if(any(y.all[,,k]>0)) TRUE else FALSE
    detected.at.all[k] <- any(y.all[, , k] > 0)
  }

  y.obs <- y.all[,,detected.at.all]    # Drop species never detected
  detected.at.site <- apply(y.obs>0, c(1,3), any)
  y.sum.all <- apply(y.all, c(1,3), sum) # Detection frequency for all species
  y.sum.obs <- y.sum.all[,detected.at.all]# Detection frequency for obs. species
  z.obs <- apply(y.all, c(1,3), max)   # Observed presence/absence matrix
  missed.sites <- z - z.obs            # Sites where species missed
  Ntotal.fs <- sum(occurring.in.sample)# Number of species in finite-sample
  Ntotal.obs <- sum(detected.at.all)   # Observed species richness (all sites)
  S.true <- apply(z, 1, sum)           # Vector of true local richness
  S.obs <- apply(z.obs, 1, sum)        # Vector of observed local richness

  # Two panels of plots
  # (1) Species-specific and community responses of occupancy to habitat
  # (2) Species-specific and community responses of detection to wind
  if(show.plot){
    par(mfrow = c(1,2), mar = c(5,5,5,3), cex.axis = 1.3, cex.lab = 1.3)
    # (1) Species-specific and community responses of occupancy to 'habitat'
    curve(plogis(beta0[1] + beta1[1] * x), -2, 2, main = "Species-specific (black) and community (red) \n response of occupancy to habitat",
    xlab = "Habitat", ylab = "Occupancy probability (psi)", ylim = c(0,1))
    for(k in 2:nspec){
      curve(plogis(beta0[k] + beta1[k] * x), -2, 2, add = TRUE)
    }
    curve(plogis(mu.lpsi + mu.beta.lpsi * x), -2, 2, col = "red", lwd = 3, add = TRUE)

    # (2) Species-specific and community responses of detection to 'wind'
    curve(plogis(alpha0[1] + alpha1[1] * x), -2, 2, main = "Species-specific (black) and community (red) \n response of detection to wind",
    xlab = "Wind", ylab = "Detection probability (p)", ylim = c(0,1))
    for(k in 2:nspec){
      curve(plogis(alpha0[k] + alpha1[k] * x), -2, 2, add = TRUE)
    }
    curve(plogis(mu.lp + mu.beta.lp * x), -2, 2, col = "red", lwd = 3, add = TRUE)

    # More plots
    # (3) True presence/absence
    # (4) Observed detection frequencies
    # (5) Sites where a species was missed
    # (6) True and observed histogram of site-specific species richness
    par(mfrow = c(2,2), cex.axis = 1.3, cex.lab = 1.3)
    mapPalette1 <- colorRampPalette(c("white", "black"))
    mapPalette2 <- colorRampPalette(c("white", "yellow", "orange", "red"))
    # (3) True presence/absence matrix (z) for all species
    # mapPalette was 2 before
    image(x = 1:nspec, y = 1:nsite, z = t(z), col = mapPalette1(4), main =
    paste("True presence/absence (z) matrix\n (finite-sample N species =",
    Ntotal.fs,")"), frame = TRUE, xlim = c(0, nspec+1),
    ylim = c(0, nsite+1), xlab = "Species", ylab = "Sites")

    # (4) Observed detection frequency for all species
    image(x = 1:nspec, y = 1:nsite, z = t(y.sum.all), col = mapPalette2(100),
      main = paste("Observed detection frequencies"),
      xlim = c(0, nspec+1), ylim = c(0, nsite+1),
      frame = TRUE, xlab = "Species", ylab = "Sites")

    # (5) Sites where a species was missed
    image(x = 1:nspec, y = 1:nsite, z = t(missed.sites), col = mapPalette1(2),
      main = paste("Matrix of missed presences\n (obs. N species =", Ntotal.obs,")"),
      frame = TRUE, xlim = c(0, nspec+1), ylim = c(0, nsite+1), xlab = "Species",
      ylab = "Sites")

    # (6) True and observed distribution of site-specific species richness
    # plot(table(S.true), col = "red", xlab = "Number of species per site",
      # xlim = c(0, max(S.true)), ylab = "Frequency",
      # main = "True (red) vs. observed (blue) \n number of species per site")
    # points(table(S.obs+(nspec/100)), col = "blue")
    histCount(S.obs, S.true,  xlab = "Number of species per site",
      main = "True (red) vs. observed (blue) \n number of species per site")
    # See file "histCount_helper.R" for details of this function.
  }

  # Output
  return(list(type=type, nsite=nsite, nrep=nrep, nspec=nspec, mean.psi=mean.psi,
    mu.lpsi=mu.lpsi, sig.lpsi=sig.lpsi, mu.beta.lpsi=mu.beta.lpsi,
    sig.beta.lpsi=sig.beta.lpsi, mean.p=mean.p, mu.lp=mu.lp, sig.lp=sig.lp,
    mu.beta.lp=mu.beta.lp, sig.beta.lp=sig.beta.lp, psi=psi, p=p, z=z,
    z.obs = z.obs, y.all=y.all, y.obs=y.obs, y.sum.all=y.sum.all,
    y.sum.obs=y.sum.obs, Ntotal.fs = Ntotal.fs, Ntotal.obs = Ntotal.obs, S.true = S.true, S.obs = S.obs))
} # endif(type=="det/nondet")

# Code for simulating community abundance data
#   (according to a community abundance model)
if(type=="counts"){
  # Prepare structures to hold data
  y.all <- y.obs <- p <- array(NA, c(nsite, nrep, nspec))
  dimnames(y.all) <- dimnames(y.obs) <- dimnames(p) <-
    list(paste("site", 1:nsite, sep=""),
    paste("rep", 1:nrep, sep=""),
    paste("sp", 1:nspec, sep=""))
  N <- lambda <- matrix(NA, nsite, nspec)
  dimnames(N) <- dimnames(lambda) <- list(paste("site", 1:nsite, sep=""),
    paste("sp", 1:nspec, sep=""))
  detected.at.all <- rep(NA, nspec)

  # Create covariates 'habitat' and 'wind'
  habitat <- sort(rnorm(nsite))     # Note 'habitat gradient' due to sorting
  wind <- matrix(rnorm(nsite * nrep), ncol=nrep)

  # Draw species-specific intercepts and slopes from their normal distributions
  # Build up linear predictors for occupancy and detection
  mu.loglam <- log(mean.lambda)
  mu.lp <- ifelse(mean.p == '1', 500, qlogis(mean.p))

  beta0 <- rnorm(nspec, mu.loglam, sig.loglam)           # lambda intercept
  beta1 <- rnorm(nspec, mu.beta.loglam, sig.beta.loglam) # lambda slope on habitat
  alpha0 <- rnorm(nspec, mu.lp, sig.lp)              # detection intercept
  alpha1 <- rnorm(nspec, mu.beta.lp, sig.beta.lp)    # detection slope on wind
  for(k in 1:nspec){
    lambda[,k] <- exp(beta0[k] + beta1[k] * habitat)
    for(j in 1:nrep){
      p[,j,k] <- plogis(alpha0[k] + alpha1[k] * wind[,j])
    }
  }

  # Distribute species over sites (simulate true abundance state)
  for(k in 1:nspec){
    N[,k] <- rpois(nsite, lambda[,k])
  }
  tmp <- apply(N, 2, sum)
  occurring.in.sample <- as.numeric(tmp > 0) # Presence/absence in study area

  # Measurement of abundance (simulate counts)
  for(k in 1:nspec) {
    for(i in 1:nsite){
      for(j in 1:nrep) {
        y.all[i,j,k] <- rbinom(1, N[i,k], p[i,j,k])
      }
    }
  detected.at.all[k] <- if(any(y.all[,,k]>0)) TRUE else FALSE
  }

  y.obs <- y.all[,,detected.at.all]    # Drop species never detected
  detected.at.site <- apply(y.obs>0, c(1,3), any)
  ymax.obs <- apply(y.all, c(1,3), max)   # Observed max count matrix
  Ntotal.fs <- sum(occurring.in.sample)# Number of species in finite-sample
  Ntotal.obs <- sum(detected.at.all)   # Observed species richness (all sites)

  # Two panels of plots
  # (1) Species-specific and community responses of lambda to habitat
  # (2) Species-specific and community responses of detection to wind
  if(show.plot){
    par(mfrow = c(1,2), mar = c(5,5,5,3), cex.axis = 1.3, cex.lab = 1.3)
    # (1) Species-specific and community responses of occupancy to 'habitat'
    curve(exp(beta0[1] + beta1[1] * x), -2, 2, main = "Species-specific (black) and community (red) \n response of lambda to habitat", xlab = "Habitat",
    ylab = "Expected abundance (lambda)")
    for(k in 1:nspec){
      curve(exp(beta0[k] + beta1[k] * x), -2, 2, add = TRUE)
    }
    curve(exp(mu.loglam + mu.beta.loglam * x), -2, 2, col = "red", lwd = 3, add = TRUE)

    # (2) Species-specific and community responses of detection to 'wind'
    curve(plogis(alpha0[1] + alpha1[1] * x), -2, 2, main = "Species-specific (black) and community (red) \n response of detection to wind",
    xlab = "Wind", ylab = "Detection probability (p)", ylim = c(0,1))
    for(k in 2:nspec){
      curve(plogis(alpha0[k] + alpha1[k] * x), -2, 2, add = TRUE)
    }
    curve(plogis(mu.lp + mu.beta.lp * x), -2, 2, col = "red", lwd = 3, add = TRUE)

    # More plots
    # (3) True abundance N (log10 + 1)
    # (4) Observed detection frequencies (log10 + 1)
    # (5) Ratio of max count to true N
    # (6) log(max count) vs. log (true N)

    par(mfrow = c(2,2), mar = c(5,6,4,2), cex.axis = 1.3, cex.lab = 1.3)
    mapPalette <- colorRampPalette(c("yellow", "orange", "red"))
    # (3) True abundance matrix (log(N+1)) for all species
    # mapPalette was 2 before
    image(x = 1:nspec, y = 1:nsite, z = log10(t(N)+1), col = mapPalette(100),
     main =   paste("True log(abundance) (log10(N)) matrix\n (finite-sample N species =", sum(occurring.in.sample),")"), frame = TRUE, xlim = c(0, nspec+1),
    zlim = c(0, log10(max(N))), xlab = "Species", ylab = "Sites")
    # (4) Observed maximum counts for all species
    image(x = 1:nspec, y = 1:nsite, z = log10(t(ymax.obs)+1), col = mapPalette(100), main = paste("Observed maximum counts (log10 + 1)"),
      xlim = c(0, nspec+1), frame = TRUE, xlab = "Species", ylab = "Sites", zlim = c(0, log10(max(N))))
    # (5) Ratio of max count to true N
    ratio <- ymax.obs/N
    ratio[ratio == "NaN"] <- 1
    image(x = 1:nspec, y = 1:nsite, z = t(ratio), col = mapPalette(100),
      main = paste("Ratio of max count to true abundance (N)"),
      xlim = c(0, nspec+1), frame = TRUE, xlab = "Species", ylab = "Sites",
      zlim = c(0, 1))

    # (6) True N and observed max count versus 'habitat'
    lims <- c(0, log10(max(N+1)))
    plot(log(N), log(ymax.obs), xlab = "True abundance (log10(N+1))",
     ylab = "Observed max count \n(log10(max+1))", xlim = lims, ylim = lims,
     main = "Observed vs. true N (log10 scale)" )
    abline(0,1)
  }

  # Output
  return(list(type=type, nsite=nsite, nrep=nrep, nspec=nspec,
    mean.lambda=mean.lambda, mu.loglam=mu.loglam, sig.loglam=sig.loglam,
    mu.beta.loglam=mu.beta.loglam, sig.beta.loglam=sig.beta.loglam,
    mean.p=mean.p, mu.lp=mu.lp, sig.lp=sig.lp, mu.beta.lp=mu.beta.lp,
    sig.beta.lp=sig.beta.lp, lambda=lambda, p=p, N=N, # y.obs=y.obs,
    y.all=y.all, y.obs=y.obs, ymax.obs=ymax.obs, Ntotal.fs=Ntotal.fs,
    Ntotal.obs=Ntotal.obs))
} # endif(type=="counts")
}

