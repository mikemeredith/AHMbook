
# AHM2 section 2.5.1, originally called 'DMsim.fn'
# A function to simulate data for a Dail-Madsen model without covariates.
# Plots added 2019-08-16, see email from Marc.

# ---------- simulator function ------------------------
simDM0 <- function(nsites = 50, nsurveys = 3, nyears = 5,
  lambda = 4, gamma = 1.5, phi = 0.8, p = 0.7, show.plots=TRUE){
  ## Simulation for multiple-visit data (from pcountOpen help file)
  ## No covariates, constant time intervals between primary periods
  # nsites: Number of sites
  # nsurveys: Number of replicate (secondary) samples within period of closure
  # nyears: Number of primary samples: years, seasons etc.
  # lambda: Initial expected abundance
  # gamma, phi: recruitment and apparent survival rates, respectively
  # p: detection probability

  # Checks and fixes for input data -----------------------------
  nsites <- round(nsites[1])
  nsurveys <- round(nsurveys[1])
  nyears <- round(nyears[1])
  stopifNegative(lambda, allowZero=FALSE)
  stopifnotProbability(phi)
  stopifnotProbability(p)
  # --------------------------------------------

  y <- array(NA, dim = c(nsites, nyears, nsurveys))
  N <- matrix(NA, nsites, nyears)
  S <- R <- matrix(NA, nsites, nyears-1)
  N[,1] <- rpois(nsites, lambda)        # Initial state
  for(t in 1:(nyears-1)) {              # State dynamics
    S[,t] <- rbinom(nsites, N[,t], phi) # Number of survivors
    R[,t] <- rpois(nsites, gamma)       # Number of recruits
    N[,t+1] <- S[,t] + R[,t]           # Number in population next year
  }
  for(j in 1:nsurveys){                   # Observation process
    y[,,j] <- rbinom(nsites*nyears, N, p)
  }
  # Put observed data into two dimensions
  yy <- array(NA, dim = c(nsites, nsurveys*nyears))
  for(t in 1:nyears){
     yy[,(nsurveys * t-(nsurveys-1)):(nsurveys*t)] <- y[,t,]
  }

  if(show.plots) {
    op <- par(mfrow = c(2,2), mar = c(5,5,4,3), cex.lab = 1.5, cex.axis = 1.5)
    on.exit(par(op))

    matplot(t(N), type = 'l', main = paste('Population trajectories under a simple DM model \nwith lambda =', lambda, ', phi =', phi, 'and gamma =', gamma, ''), lty = 1, lwd = 3, las = 1, frame = FALSE, xlab = 'Year', ylab = 'N')
    matplot(t(S), type = 'l', main = 'Number of apparent survivors', lty = 1, lwd = 3, las = 1, frame = FALSE, xlab = 'Year', ylab = 'S')
    hist(N[,1], main = 'Distribution of N in first year', breaks = 50, col = 'grey')
    hist(N[,nyears], main = 'Distribution of N in last year', breaks = 50, col = 'grey')
  }
  return(list(
    # -------------- arguments input -------------------
    nsites = nsites, nsurveys = nsurveys, nyears = nyears, lambda = lambda,
    gamma = gamma, phi = phi, p = p,
    # ----------- values generated -------------------------
    N = N,          # true number of individuals, nsites x nyears
    S = S, R = R,   # number of survivors, recruits, nsites x (nyears-1)
    y = y,          # number detected, nsites x nyears x nsurveys
    yy = yy))       # number detected as a 2D matrix, nsites x (nyears*nsurveys)
}
# -------------------- end function ----------------------------
