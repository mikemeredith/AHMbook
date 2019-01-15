
# 13.5.5 A function to simulate data for a Dail-Madsen model with covariates.

simDM <- function(nsite = 50, nsurvey = 3, nyear = 5,
  mean.lambda = 4, mean.gamma.rel = 0.5,
  mean.phi = 0.8, mean.p = 0.7,
  beta.lam = 1, beta.gamma = 1, beta.phi = -1, beta.p = -1,
  show.plot=TRUE){
  # Simulation for multiple-visit data,
  # constant time intervals between primary periods
  # nsite: Number of sites
  # nsurvey: Number of replicate (secondary) samples within period of closure
  # nyear: Number of primary samples: years, seasons etc.
  # mean.lambda: Initial expected abundance at cov.lam = 0
  # mean.gamma.rel, mean.phi: recruitment and apparent survival rates,
  #   respectively, at values of cov.gamma and cov.phi equal to 0
  # mean.p: detection probability at cov. p = 0
  # beta.X is the slope of parameter X (link transformed) on the respective covariate
  y <- p <- array(NA, dim = c(nsite, nyear, nsurvey))
  N <- matrix(NA, nsite, nyear)
  S <- R <- matrix(NA, nsite, nyear-1)
  cov.lam <- runif(nsite, -1, 1)
  cov.gamma <- runif(nsite, -1, 1)
  cov.phi <- runif(nsite, -1, 1)
  cov.p <- array(runif(nsite*nyear*nsurvey, -1, 1), dim = dim(y))

  lambda <- exp(log(mean.lambda) + beta.lam * cov.lam)
  N[,1] <- rpois(nsite, lambda)        # Initial state

  phi <- plogis(qlogis(mean.phi) + beta.phi * cov.phi)
  gamma <- exp(log(mean.gamma.rel) + beta.gamma * cov.gamma)

  for(t in 1:(nyear-1)) {              # State dynamics
     S[,t] <- rbinom(nsite, N[,t], phi)
     R[,t] <- rpois(nsite, N[,(t)]*gamma)    # Simulate in 'relative' mode
     N[,t+1] <- S[,t] + R[,t]
  }
  for(i in 1:nsite){                   # Observation process
     for(t in 1:nyear){
        for(j in 1:nsurvey){
           p[i,t,j] <- plogis(qlogis(mean.p) + beta.p * cov.p[i,t,j])
           y[i,t,j] <- rbinom(1, N[i,t], p[i,t,j])
        }
     }
  }

  # Put observed data into two dimensions
  yy <- ccov.p <- array(NA, dim = c(nsite, nsurvey*nyear))
  for(t in 1:nyear){
    yy[,(nsurvey * t-(nsurvey-1)):(nsurvey*t)] <- y[,t,]
    ccov.p[,(nsurvey * t-(nsurvey-1)):(nsurvey*t)] <- cov.p[,t,]
  }

  # Visualisations
  if(show.plot) {
    op <- par(mfrow = c(3, 2), mar = c(5,5,4,3), cex.lab = 1.5, cex.axis = 1.5)
    on.exit(par(op))

    matplot(t(N), type = 'l',
      main = paste('Population trajectories under a simple DM model \nwith mean lambda =',
      mean.lambda, ', mean gamma =', mean.gamma.rel, ' and mean phi =', mean.phi, ''),
      lty = 1, lwd = 3, las = 1, frame = FALSE, xlab = 'Year', ylab = 'N')
    matplot(t(S), type = 'l', main = 'Number of apparent survivors', lty = 1, lwd = 3, las = 1,
      frame = FALSE, xlab = 'Year', ylab = 'Survivors (S)')
    matplot(t(R), type = 'l', main = 'Number of recruits', lty = 1, lwd = 3, las = 1,
      frame = FALSE, xlab = 'Year', ylab = 'Recruits (R)')
    matplot(t(apply(p, c(1,2), mean)), type = 'l',
      main = 'Average detection probability per site and year', lty = 1, lwd = 3, las = 1,
      frame = FALSE, xlab = 'Year', ylab = 'Average p')
    hist(N[,1], main = 'Distribution of N in first year', breaks = 50, col = 'grey')
    hist(N[,nyear], main = 'Distribution of N in last year', breaks = 50, col = 'grey')
  }

  # Output
  return(list(
    # -------------- arguments input -------------------
    nsite = nsite, nsurvey = nsurvey, nyear = nyear, mean.lambda = mean.lambda, mean.gamma.rel = mean.gamma.rel, mean.phi = mean.phi, mean.p = mean.p, beta.lam = beta.lam, beta.gamma = beta.gamma, beta.phi = beta.phi, beta.p = beta.p,
    # ----------- values generated -------------------------
    cov.lam = cov.lam, cov.gamma = cov.gamma, cov.phi = cov.phi,  # covariates
    cov.p = cov.p,   # covariate for p, nsite x nyear x nsurvey
    ccov.p = ccov.p, # covariate for p as 2D matrix, nsite x (nyear*nsurvey)
    N = N,           # true number of individuals, nsite x nyear
    S = S, R = R,    # number of survivors, recruits, nsite x (nyear-1)
    p = p,           # probability of detection, nsite x nyear x nsurvey
    y = y,           # number detected, nsite x nyear x nsurvey
    yy = yy))        # number detected as a 2D matrix, nsite x (nyear*nsurvey)
}


