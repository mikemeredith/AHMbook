
# AHM2 section 2.7.1 Simulate data for a multinomial-mixture model

# Called data.fn in the draft
# with R=nsites, J=nsurveys, K=nyears

simMultMix <- function(nsites = 100, nsurveys = 3, nyears = 4,
  lambda = 3, theta = 0.5, p = 0.3){
  # Simulate data using the multinomial-Poisson model with a
  # repeated constant-interval removal design (written by R.B. Chandler)
  # R :        Number of sites (variable)
  # K <- 4     # Number of primary periods
  # J <- 3     # Number of secondary periods
  # lambda, theta and p: expected abundance, availability and detection prob.

  # Checks and fixes for input data -----------------------------
  nsites <- round(nsites[1])
  nsurveys <- round(nsurveys[1])
  nyears <- round(nyears[1])
  stopifNegative(lambda, allowZero=FALSE)
  stopifnotProbability(theta)
  stopifnotProbability(p)
  # --------------------------------------------

  y <- array(NA, c(nsites, nyears, nsurveys))
  M <- rpois(nsites, lambda)   # Local population size
  N <- matrix(NA, nsites, nyears)   # Individuals available for detection

  for(i in 1:nsites) {
     N[i,] <- rbinom(nyears, M[i], theta)
     y[i,,1] <- rbinom(nyears, N[i,], p)    # Observe some
     Nleft1 <- N[i,] - y[i,,1]     # Remove them
     y[i,,2] <- rbinom(nyears, Nleft1, p)   # ...
     Nleft2 <- Nleft1 - y[i,,2]
     y[i,,3] <- rbinom(nyears, Nleft2, p)
  }
  y2d <- cbind(y[,1,], y[,2,], y[,3,], y[,4,])
  return(list(
      # ------ input arguments ------
      nsites = nsites, nsurveys = nsurveys, nyears = nyears, lambda = lambda,
      theta = theta, p = p,
      # ------ generated values ------
      M = M, # local population size
      N = N, # animals available for detection
      y = y, # sites x years x surveys array of observations
      y2d = y2d)) # y in 2d format, sites x (years*surveys) matrix
}
