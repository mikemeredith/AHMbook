
# Function for section 2.9.1, Diseased Frogs

simFrogDisease <- function(nsites = 100, nyears = 3, nsurveys = 3,
                    alpha.lam = 3,             # Mean abundance at t=1
                    omega = c(0.9, 0.7),       # State-specific survival
                    gamma = c(2,1),            # State-specific recruitment
                    p = c(0.8, 0.8, 0.8),      # Detection probability
                    recovery = 0.1,        # Pr recovery given diseased
                    infection = 0.1){      # Pr infection given not diseased

  # Empty matrices to hold the data
  yN <- yI <- array(NA, dim = c(nsites, nyears, nsurveys))
  NI <- NN <- array(NA, dim = c(nsites, nsurveys))

  # First season
  NN[,1] <- rpois(n = nsites, lambda = alpha.lam)
  NI[,1] <- rpois(n = nsites, lambda = alpha.lam)
  for(i in 1:nsites){
      for(j in 1:nyears){
        yN[i,j, 1] <- rbinom(n = 1, NN[i,1], p[1])
        yI[i,j, 1] <- rbinom(n = 1, NI[i,1], p[1])
      }
  }
  SN <- SI <- GI <- GN <- TrN <- TrI <- array(0, dim = c(nsites, nsurveys-1))

  # Second and subsequent seasons
  for(k in 2:nsurveys){
    for(i in 1:nsites){
      if(NN[i,k-1]>0){
        SN[i, k-1] <- rbinom(n=1, size=NN[i,k-1], prob=omega[1])	# Survival of uninfecteds
        TrN[i,k-1] <- rbinom(n=1, size=SN[i,k-1], prob=infection)   # Getting infected - lost from NN, and gained by NI
      }
      if(NI[i,k-1]>0){
        SI[i, k-1] <-  rbinom(n=1, size=NI[i,k-1], prob=omega[2])   # Survival of infecteds
        TrI[i, k-1] <- rbinom(n=1, size=SI[i,k-1], prob=recovery) 	# Losing infection - lost from NI and gained by NN
      }
      # Recruitment
      GI[i, k-1] <- rpois(1, lambda = gamma[2])
      GN[i, k-1] <- rpois(1, lambda = gamma[1])
    }
    # Total population size
    NI[,k] <-  SI[,k-1] + GI[,k-1]  + TrN[,k-1] - TrI[,k-1]
    NN[,k] <-  SN[,k-1] + GN[,k-1]  + TrI[,k-1] - TrN[,k-1]
  }
  for(i in 1:nsites){
    for(j in 1:nyears){
      for(k in 2:nsurveys){
            yN[i, j, k] <- rbinom(n = 1, NN[i,k], p[k])
            yI[i, j, k] <- rbinom(n = 1, NI[i,k], p[k])
      }
    }
  }
  return(list(
    # --------------- arguments input ------------------------
    nsites = nsites, nyears = nyears, nsurveys = nsurveys,alpha.lam= alpha.lam,omega = omega,gamma = gamma,
    infection = infection, recovery = recovery,
    # ---------------- generated values ----------------------
    SN = SN,   # sites x intervals, number of noninfected frogs surviving
    SI = SI,   # sites x intervals, number of infected frogs surviving
    GN = GN,   # sites x intervals, number of noninfected frogs recruited
    GI = GI,   # sites x intervals, number of infected frogs recruited
    TrI = TrI, # sites x intervals, number of infected frogs recovering
    TrN = TrN, # sites x intervals, number of noninfected frogs becoming infected
    NN = NN,   # sites x years, number of noninfected frogs in the population
    NI = NI,   # sites x years, number of infected frogs in the population
    p = p,     # length nyears, probability of detection
    yN = yN,   # sites x years x surveys, number of noninfected frogs detected
    yI = yI))  # sites x years x surveys, number of infected frogs detected
}
