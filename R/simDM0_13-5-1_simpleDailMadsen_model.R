
# section 13.5.1, originally called 'DMsim.fn'
# A function to simulate data for a Dail-Madsen model without covariates.

# ---------- simulator function ------------------------
simDM0 <- function(nsite = 50, nsurvey = 3, nyear = 5, 
  lambda = 4, gamma = 1.5, phi = 0.8, p = 0.7){
  ## Simulation for multiple-visit data (from pcountOpen help file)
  ## No covariates, constant time intervals between primary periods
  # nsite: Number of sites
  # nsurvey: Number of replicate (secondary) samples within period of closure
  # nyear: Number of primary samples: years, seasons etc.
  # lambda: Initial expected abundance
  # gamma, phi: recruitment and apparent survival rates, respectively
  # p: detection probability

  y <- array(NA, dim = c(nsite, nyear, nsurvey))
  N <- matrix(NA, nsite, nyear)
  S <- R <- matrix(NA, nsite, nyear-1)
  N[,1] <- rpois(nsite, lambda)        # Initial state
  for(t in 1:(nyear-1)) {              # State dynamics
    S[,t] <- rbinom(nsite, N[,t], phi) # Number of survivors
    R[,t] <- rpois(nsite, gamma)       # Number of recruits
    N[,t+1] <- S[,t] + R[,t]           # Number in population next year
  }
  for(j in 1:nsurvey){                   # Observation process
    y[,,j] <- rbinom(nsite*nyear, N, p)
  }
  # Put observed data into two dimensions
  yy <- array(NA, dim = c(nsite, nsurvey*nyear))
  for(t in 1:nyear){
     yy[,(nsurvey * t-(nsurvey-1)):(nsurvey*t)] <- y[,t,]
  }
  return(list(
    # -------------- arguments input -------------------
    nsite = nsite, nsurvey = nsurvey, nyear = nyear, lambda = lambda,
    gamma = gamma, phi = phi, p = p, 
    # ----------- values generated -------------------------
    N = N,          # true number of individuals, nsites x nyears
    S = S, R = R,   # number of survivors, recruits, nsites x (nyears-1)
    y = y,          # number detected, nsites x nyears x nsurvey
    yy = yy))       # number detected as a 2D matrix, nsites x (nyears*nsurvey)
}
# -------------------- end function ----------------------------
