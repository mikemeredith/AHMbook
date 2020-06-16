

# AHM2 section 1.7.1 Simulation of a demographic state-space model

# Called simpop in the draft

# Define function to simulate the data
# ----------------- Start function definition -------------------
simPOP <- function(
  M = 100,                      # number of sites
  T = 10,                       # number of years
  mean.lam = 3,                 # mean abundance for year 1
  beta.lam = 0,                 # covariate coefficient for lambda
  sd.log.lam = 0,               # over-dispersion in lambda
  mean.gamma = 1.0,             # mean population growth rate
  beta.gamma = 0,               # covariate coefficient for gamma
  sd.log.gamma.site = 0,        # SD of site effects
  sd.log.gamma.time = 0,        # SD of time effects
  sd.log.gamma.survey = 0,      # SD of survey (site+time) effects
  sd.rho = 0,                   # random immigration term
  mean.p = 0.6,                 # mean detection probability
  beta.p = 0,                   # covariate coefficient for p
  sd.logit.p.site = 0,          # SD of site effects
  sd.logit.p.time = 0,          # SD of time effects
  sd.logit.p.survey = 0,        # SD of survey effects
  show.plot = TRUE){         # controls plotting

  # Simulate multiple time-series of counts under a pure Markov model (with exponential population model) or
  # under an extended Markov model (with exponential-plus-random-immigration population model;
  # see Sollmann et al., Ecology, 2015)
  # Default is Markov model, setting sd.rho to a value greater than 0 changes to extended Markov and sets the amount of random immigration.

  # Checks and fixes for input data -----------------
  M <- round(M[1])
  T <- round(T[1])
  stopifNegative(mean.lam)
  stopifNegative(sd.log.lam)
  stopifNegative(sd.log.gamma.site)
  stopifNegative(sd.log.gamma.time)
  stopifNegative(sd.log.gamma.survey)
  stopifNegative(sd.rho)
  stopifnotProbability(mean.p)
  stopifNegative(sd.logit.p.site)
  stopifNegative(sd.logit.p.time)
  stopifNegative(sd.logit.p.survey)
  # -----------------------------------------------------

  # Create arrays needed (for observed counts, latent states, gamma, p
  # C <- N <- gamma <- p <- array(NA, dim = c(M, T))
  # (Could also have this: gamma <- array(NA, dim = c(M, T-1))) Mike prefers this!
  C <- N <- p <- array(NA, dim = c(M, T))
  gamma <- array(NA, dim = c(M, T-1))

  # Assemble lambda
  Xsite1 <- runif(M, -1, 1)                # Site covariate that affects initial abundance
  eps.N <- rnorm(M, 0, sd.log.lam)         # Site overdispersion at t = 1
  lambda <- exp(log(mean.lam) + beta.lam * Xsite1 + eps.N)

  # Assemble gamma (last column will remain NA)
  Xsiteyear1 <- matrix(runif(M*T, -1, 1), nrow = M, ncol = T) # Yearly site covariate that affects gamma
  eps.gamma.site <- rnorm(M, 0, sd.log.gamma.site)      # spatial (= site) effects in gamma
  eps.gamma.time <- rnorm(T, 0, sd.log.gamma.time)      # temporal (=primary occasion) effects in gamma
  eps.gamma.survey <- matrix(rnorm(M*T, 0, sd.log.gamma.survey), nrow = M, ncol = T) # survey effects in gamma (i.e., site by occasion)
  for(t in 1:(T-1)){
    gamma[,t] <- exp(log(mean.gamma) + beta.gamma * Xsiteyear1[,t] + eps.gamma.site + eps.gamma.time[t] + eps.gamma.survey[,t])
  }

  # Draw values of random immigration rate (rho)
  if(sd.rho == 0){  # Markovian dynamics
    rho <- rep(0, T-1)
  }
  if(sd.rho > 0){   # Extended Markovian dynamics
    logrho <- rnorm(T-1, 0, sd.rho)
    rho <- exp(logrho)
  }

  # Assemble p
  Xsiteyear2 <- matrix(runif(M*T, -1, 1), nrow = M, ncol = T) # Yearly site covariate that affects p
  eps.p.site <- rnorm(M, 0, sd.logit.p.site)             # spatial (= site) effects in p
  eps.p.time <- rnorm(T, 0, sd.logit.p.time)           # temporal (=primary occasion) effects in p
  eps.p.survey <- matrix(rnorm(M*T, 0, sd.logit.p.survey), nrow = M, ncol = T) # survey effects in p (i.e., site by occasion)
  for(t in 1:T){
    p[,t] <- plogis(qlogis(mean.p) + beta.p * Xsiteyear2[,t] + eps.p.site + eps.p.time[t] + eps.p.survey[,t])
  }

  # Simulate initial state
  N[,1] <- rpois(M, lambda)

  # Simulate later states
  for(t in 2:T){
    N[,t] <- rpois(M, N[,t-1] * gamma[,t-1] + rho[t-1])
  }

  # Simulate binomial observation
  for(t in 1:T){
    C[,t] <- rbinom(M, N[,t], p[,t])
  }

  # Tally up number of extinct populations and compute extinction rate
  Nextinct <- sum(N[,T] == 0)
  extrate <- Nextinct / M

  # Tally up number of years in which a pop is at zero
  zeroNyears <- sum(N == 0)

  # Add up annual total population size across all sites
  sumN <- apply(N, 2, sum)

  # Compute realized population growth rate based on total, realized population size per year
  gammaX <- numeric(T-1)
  for(t in 2:T){
    gammaX[t-1] <- sumN[t] / sumN[t-1]
  }

  # Graphical output
  if(show.plot) {
    # Restore graphical settings on exit
    oldpar <- par("mfrow", "mar")
    oldAsk <- devAskNewPage(ask = dev.interactive(orNone=TRUE))
    on.exit({par(oldpar); devAskNewPage(oldAsk)})

    tryPlot <- try( {
      par(mfrow = c(1,3))
      hist(lambda, breaks = 100, main = 'lambda', col = 'grey')
      hist(gamma, breaks = 100, main = 'gamma', col = 'grey')
      hist(p, breaks = 100, main = 'p', col = 'grey')

      par(mfrow = c(1, 3))
      hist(N, breaks = 100, main = 'N', col = 'grey')
      hist(C, breaks = 100, main = 'C', col = 'grey')
      plot(N, C, xlab = 'True N', ylab = 'Observed C', frame = FALSE)
      abline(0,1)

      par(mfrow = c(2, 2))
      ylim <- range(c(N, C))
      matplot(t(N), type = 'l', lty = 1, main = 'Trajectories of true N', frame = FALSE, ylim = ylim)
      matplot(t(C), type = 'l', lty = 1, main = 'Trajectories of observed C', frame = FALSE, ylim = ylim)
      plot(table(N[,1]), main = 'Initial N', frame = FALSE)
      plot(table(N[,T]), main = 'Final N', frame = FALSE)
    }, silent = TRUE)
    if(inherits(tryPlot, "try-error"))
      tryPlotError(tryPlot)
  }
  # Numeric output
  return(list(
    # ----------------- arguments input ----------------------
    M = M, T = T, mean.lam = mean.lam, beta.lam = beta.lam, sd.log.lam = sd.log.lam, mean.gamma = mean.gamma, beta.gamma = beta.gamma, sd.log.gamma.site = sd.log.gamma.site, sd.log.gamma.time = sd.log.gamma.time, sd.log.gamma.survey = sd.log.gamma.survey, sd.rho = sd.rho, mean.p = mean.p, beta.p = beta.p, sd.logit.p.site = sd.logit.p.site, sd.logit.p.time = sd.logit.p.time, sd.logit.p.survey = sd.logit.p.survey,
    # ------------------ generated values ----------------------
    Xsite1 = Xsite1,          # M vector, site covariate affecting initial abundance
    Xsiteyear1 = Xsiteyear1,  # MxT matrix, yearly site covariate affecting gamma
    Xsiteyear2 = Xsiteyear2,  # MxT matrix, yearly site covariate affecting p
    eps.N = eps.N,            # M vector, site overdispersion at t = 1
    lambda = lambda,          # M vector, abundance in year 1
    eps.gamma.site = eps.gamma.site,     # M vector, random site effect for gamma
    eps.gamma.time = eps.gamma.time,     # T vector, random time effect for gamma
    eps.gamma.survey = eps.gamma.survey, # MxT matrix, random survey effect for gamma
    gamma = gamma,                       # MxT matrix, population growth rate
    rho = rho,                           # T-1 vector, immigration rate
    eps.p.site = eps.p.site,     # M vector, random site effect for p
    eps.p.time = eps.p.time,     # T vector, random time effect for p
    eps.p.survey = eps.p.survey, # MxT matrix, random survey effect for p
    p = p,                       # MxT matrix, detection probability
    N = N,                       # MxT matrix, true population
    C = C,                       # MxT matrix, counts
    zeroNyears = zeroNyears,     # scalar, sum(N == 0)
    Nextinct = Nextinct,         # scalar, number of sites where N ==0 at time T
    extrate = extrate,           # scalar, proportion of sites where N ==0 at time T
    sumN = sumN,                 # T vector, total population in each year
    gammaX = gammaX))            # T-1 vector, realized population growth rate
}
# ----------------- End function definition -------------------
