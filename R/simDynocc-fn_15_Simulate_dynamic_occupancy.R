# Functions for the book Applied Hierarchical Modeling in Ecology (AHM)
# Marc Kéry & Andy Royle, Academic Press, 2016.

# simDynocc.fn - section 15.??


# Function to generate data under a non-spatial dynamic occupancy model

simDynocc.fn <- function(M = 250, J = 3, K = 10, mean.psi1 = 0.4, beta.Xpsi1 = 0,
range.phi = c(0.5, 1), beta.Xphi = 0, range.gamma = c(0, 0.5), beta.Xgamma = 0, range.p = c(0.1, 0.9), beta.Xp = 0, range.sd.site = c(0, 0), range.sd.survey = c(0, 0), seed = Sys.time()) {
#
# Written by Marc Kéry, 4 Dec 2014
#
# Function to simulate detection/nondetection data under a very general
#     dynamic site-occ model, including:
#   - annual variation in the probabilities of patch persistence, colonization
#     and detection is specified by the bounds of a uniform distribution.
#   - one site-specific (???) covariate each is allowed to affect
#     every parameter
#   - autologistic effects in persistence and colonization probability can be
#     chosen, which fits a logistic regression of these parameters on the
#     proportion of occupied neighbouring cells (in a queen's or
#     2nd order neighbourhood)
#   - Additional detection heterogeneity at the site- or the survey level,
#     with the possibility of a temporal trend in this heterogeneity. E.g.,
#     an annual trend in detection heterogeneity at the site or
#     the survey level is specified by the value in the first and the last year.
#     Hence, range.sd.site = c(0, 1) will result in a linear trend in the
#     magnitude of site heterogeneity in detection from 0 in the first year to
#     1 in the last year.

# Function arguments:
# -------------------
# LATER: need number pixels per side of square ------------------
# M – Number of sites
# J – Number of replicate surveys within a year (= season)
# K – Number of years (or 'seasons')
# mean.psi1 – average occupancy probability in first year
# range.p – bounds of uniform distribution from which annual p drawn
# range.psi and range.gamma – same for survival and colonization probability
# -------------------
# beta.Xpsi1, beta.Xphi, beta.Xgamma, beta.Xp – coefficients of
# environmental covariates in probabilities of initial occupancy, persistence,
# colonization and detection.
# spatial corr : LATER: need those for each covariate field --------
# -------------------
# range.sd.site: sd of normal distribution to model logit-normal noise in p
# at the site level in the first and the last year of the simulation
# range.sd.survey: sd of normal distribution to model logit-normal noise in p
# at the site/year/rep = ‘survey’ level, in the first and the last year
#
# For the sd and error.rate arguments, if the two values in the range are the
# same, a constant value is assumed over time, while if they are different, a
# linear trend is assumed over time.

set.seed(as.numeric(seed))


# Set up arrays needed
site <- 1:M                             # Sites
year <- 1:K                             # Years
psi <- muZ <- z <- array(dim = c(M, K)) # Occupancy, occurrence
phi <- gamma <- array(NA, dim = c(M, (K-1))) # Survival, colonisation
y <- p <- array(NA, dim = c(M, J, K))        # Det. histories and p


# Create covariates
Xpsi1 <- runif(M, -2, 2)                     # Site covariate for psi1
Xphi <- Xgamma <- array(runif(M*K, -2, 2), dim = c(M,K))
                                             # Yearly-site covariates
Xp <- array(runif(M*J*K,-2,2),dim=c(M,J,K))  # Observ. cov.

# (1) Simulate all parameter values
# (a) State process parameters
psi[,1] <- plogis(logit(mean.psi1) + beta.Xpsi1 * Xpsi1) # psi1
mean.phi <- runif(n = K-1, min = range.phi[1], max = range.phi[2])
mean.gamma <- runif(n = K-1, min = range.gamma[1], max = range.gamma[2])
for(k in 1:(K-1)){
   phi[,k] <- plogis(logit(mean.phi[k]) + beta.Xphi * Xphi[,k])
   gamma[,k] <- plogis(logit(mean.gamma[k]) + beta.Xgamma * Xgamma[,k])
}

# (b) Observation process parameters
mean.p <- runif(n = K, min = range.p[1], max = range.p[2])
sd.site <- seq(from = range.sd.site[1], to = range.sd.site[2], length.out = K)
sd.survey <- seq(from = range.sd.survey[1], to = range.sd.survey[2], length.out = K)
for(i in 1:M){
  for(k in 1:K){
    eps1 <- rnorm(n = M, mean = 0, sd = sd.site[k])   # Site random eff.
    eps2 <- rnorm(n = J, mean = 0, sd = sd.survey[k]) # Survey random eff.
    for(j in 1:J){
      p[i,j,k] <- plogis(logit(mean.p[k])+beta.Xp*Xp[i,j,k]+eps1[i]+eps2[j])
    }
  }
}


# (2) Simulate the true system dynamics (state process)
# First year
z[,1] <- rbinom(M, 1, psi[,1])   # Initial occurrence state
# Years 2:K
for(i in 1:M){                   # Loop over sites
   for(k in 2:K){                # Loop over years
      muZ[i,k] <- z[i, k-1]*phi[i,k-1] + (1-z[i, k-1])*gamma[i,k-1]
      z[i,k] <- rbinom(1, 1, muZ[i,k])
   }
}

# (3) Simulate observation process
for(i in 1:M){                     # Loop over sites
   for(k in 1:K){                  # Loop over years
      for(j in 1:J){               # Loop over replicates
         prob <- z[i,k] * p[i,j,k] # zero out p for unoccupied sites
         y[i,j,k] <- rbinom(1, 1, prob)
      }
   }
}

# (4) Plots and computation of derived quantities
# Plot realised occupancy
plot(year, apply(z, 2, mean), type = "l", xlab = "Year", ylab = "Occupancy or Detection prob.", col = "red", xlim = c(0,K+1), ylim = c(0,1), lwd = 2, lty = 1, frame.plot = FALSE, las = 1)
lines(year, mean.p , type = "l", col = "red", lwd = 2, lty = 2)
# Plot apparent occupancy
psi.app <- apply(apply(y, c(1,3), max), 2, mean)
lines(year, psi.app, type = "l", col = "black", lwd = 2)
text(0.85*K, 0.1, labels = "red solid – true occupancy\n red dashed – detection\n black – observed occupancy", cex = 0.7)


# Compute annual population occupancy
for(i in 1:M){
   for (k in 2:K){
      psi[i,k] <- psi[i,k-1]*phi[i,k-1] + (1-psi[i,k-1])*gamma[i,k-1]
   }
}
mean.psi <- apply(psi, 2, mean)              # Average psi over sites

# Return data
return(list(M=M, J=J, K=K, mean.psi1=mean.psi1, beta.Xpsi1=beta.Xpsi1,
range.phi=range.phi, beta.Xphi=beta.Xphi, range.gamma=range.gamma, beta.Xgamma=beta.Xgamma, range.p=range.p, beta.Xp=beta.Xp, range.sd.site=range.sd.site, range.sd.survey=range.sd.survey, sd.site=sd.site, sd.survey=sd.survey, mean.phi=mean.phi, mean.gamma=mean.gamma, mean.p=mean.p, psi=psi, mean.psi=mean.psi, psi.app=psi.app, z=z, phi=phi, gamma=gamma, p=p, y = y, seed=as.numeric(seed)))
}
