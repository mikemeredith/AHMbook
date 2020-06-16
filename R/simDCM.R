
# 16.2 A general function to simulate data under the DCM model

simDCM <- function(nspecies = 50, nsites = 100, nsurveys = 3, nyears = 10,
  mean.psi1 = 0.4, sig.lpsi1 = 1, mu.beta.lpsi1 = 0, sig.beta.lpsi1 = 0,
  range.mean.phi = c(0.8, 0.8), sig.lphi = 1,
  mu.beta.lphi = 0, sig.beta.lphi = 0,
  range.mean.gamma = c(0.2, 0.2), sig.lgamma = 1,
  mu.beta.lgamma = 0, sig.beta.lgamma = 0,
  range.mean.p = c(0.5, 0.5), sig.lp = 1,
  mu.beta.lp = 0, sig.beta.lp = 0,
  range.beta1.survey = c(0, 0), range.beta2.survey = c(0, 0),
  trend.sd.site = c(0, 0), trend.sd.survey = c(0, 0),
  show.plot = TRUE, verbose = TRUE) {
#
# Written by Marc Kery, 28 Nov 2016
#
# Function is based on the dynocc function 'simDynocc' (AHM2, Chap 15) and
# on the community occupancy function 'simComm' (AHM1, Chap 11).
#
# Function to simulate detection/nondetection data under a general
#     dynamic community (site-occ) model, including:
#   * annual variation in the probabilities of patch persistence, colonization
#     and detection is specified by the bounds of a uniform distribution.
#   * species heterogeneity around the means is specified by the SD of a normal
#     distribution and expressed on the logit scale
#   * one covariate is allowed to a parameter (site covariate for psi1,
#     site-year covariate for phi and gamma and site-year-rep for p).
#     Each covariate is allowed to differ among species again according
#     to a logit-normal model of heterogeneity.
#   * Additional detection heterogeneity at the site- or the survey level,
#     with the possibility of a temporal trend in this heterogeneity. E.g.,
#     an annual trend in detection heterogeneity at the site or
#     the survey level is specified by the value in the first and the last year.
#     Hence, range.sd.site = c(0, 1) will result in a linear trend in the
#     magnitude of site heterogeneity in detection from 0 in the first year to
#     1 in the last year.
#   * Additional detection heterogeneity that varies over the survey (= occasion)
#     according to a quadratic effect of occasion number (to model phenology of
#     an insect species for instance).
#   * These last two types of detection heterogeneity are not (yet) allowed
#     to be species-specific.
#
# Function arguments:
# -------------------
# *** Sample size arguments ***
# nspecies - Number of species (typically called N in AHM book)
# nsites - Number of sites (M)
# nsurveys - Number of replicate surveys within a year (= season) (J)
# nyears - Number of years (or 'seasons') (T)
#
# *** Arguments for mean parameters for the intercepts ***
# mean.psi1 - average occupancy probability in first year
# range.mean.p - bounds of uniform distribution from which annual p drawn
# range.mean.phi and range.mean.gamma - same for persistence and colonization prob.
# -------------------
# *** Arguments for mean parameters for the slopes ***
# mu.beta.lpsi1, mu.beta.lphi, mu.beta.lgamma, mu.beta.lp - coefficients of
#      covariates in probabilities of initial occupancy, persistence,
#      colonization and detection. These are the probability-scale means of
#      the normal distibutions, from which species-specific slopes are drawn
# -------------------
# *** Args. for species-specific heterogeneity in intercepts and slopes ***
# sig.lpsi1: sd of the normal distribution from which species-specific occupancy
#      intercepts are drawn (centered on logit(mean.psi1)), on logit scale
# sig.beta.lpsi1: sd of the normal distribution from which species-specific
#      slopes are drawn (centered on mu.beta.lpsi1)
# sig.lphi: sd of the normal distribution from which species-specific persistence
#      intercepts are drawn (centered on logit(mean.phi), which are year-specific),
#      on logit scale
# sig.beta.lphi: sd of the normal distribution from which species-specific
#      persistence slopes are drawn (centered on mu.beta.lphi)
# sig.lgamma: sd of the normal distribution from which species-specific
#      colonization intercepts are drawn (centered on logit(mean.gamma),
#      which are year-specific), on logit scale
# sig.beta.lgamma: sd of the normal distribution from which species-specific
#      colonization slopes are drawn (centered on mu.beta.lgamma)
# sig.lp: sd of the normal distribution from which species-specific
#      detection intercepts are drawn (centered on logit(mean.p),
#      which are year-specific), on logit scale
# sig.beta.lp: sd of the normal distribution from which species-specific
#      detection slopes are drawn (centered on mu.beta.lp)
# -------------------
# *** Args. for detection heterogeneity among sites and surveys
#      (this part of the model is NOT species-specific) ***
# trend.sd.site: sd of normal distribution to model logit-normal noise in p
#      at the site level in the first and the last year of the simulation.
# trend.sd.survey: sd of normal distribution to model logit-normal noise in p
#      at the site/year/rep = 'survey' level, in the first and the last year
# For the sd and error.rate arguments, if the two values in the range are the
#      same, a constant value is assumed over time, while if they are different,
#      a linear trend is assumed over time.
# -------------------
# *** Args. for detection heterogeneity among occasions within a season
#      (this part of model again NOT species-specific)
# range.beta1.survey is the range of the annual variation in the linear effect
#     of survey (i.e., of month 1-12) on the product of
#     availability and detection linear and quadratic effect of survey
# range.beta2.survey is the same for the quadratic effect of survey
#
# show.plot: if TRUE, plots are produced. Usually set to FALSE when running sims.
#

# Checks and fixes for input data -----------------------------
nspecies <- round(nspecies[1])
nsites <- round(nsites[1])
nsurveys <- round(nsurveys[1])
nyears <- round(nyears[1])
stopifnotGreaterthan (nyears, 1)
stopifnotProbability(mean.psi1)
stopifNegative(sig.lpsi1, allowZero=TRUE)
# mu.beta.lpsi1
stopifNegative(sig.beta.lpsi1, allowZero=TRUE)
stopifnotProbability(range.mean.phi) # bounds
stopifNegative(sig.lphi, allowZero=TRUE)
stopifNegative(sig.beta.lphi, allowZero=TRUE)
stopifnotProbability(range.mean.gamma) # bounds
stopifNegative(sig.lgamma, allowZero=TRUE)
stopifNegative(sig.beta.lgamma, allowZero=TRUE)
stopifnotProbability(range.mean.p) # bounds
stopifNegative(sig.lp, allowZero=TRUE)
stopifNegative(sig.beta.lp, allowZero=TRUE)
# --------------------------------------------

# Set up arrays needed
spec <- 1:nspecies                         # Species
site <- 1:nsites                           # Sites
year <- 1:nyears                           # Years
# visit <- 1:nsurveys                      # Visit
# month <- 1:nsurveys                           # Months (= surveys)
survey <- 1:nsurveys                           # Months (= surveys)
psi <- muZ <- z <- array(dim = c(nsites, nyears, nspecies), dimnames =
   list(paste('Site', site, sep = ''), paste('Year', year, sep = ''),
   paste('Spec', spec, sep = ''))) # Occupancy, occurrence
phi <- gamma <- array(NA, dim = c(nsites, (nyears-1), nspecies), dimnames =
   list(paste('Site', site, sep = ''), paste('Year', year[-nyears], sep = ''),
   paste('Spec', spec, sep = ''))) # Survival, colonisation
y <- p <- array(NA, dim = c(nsites, nsurveys, nyears, nspecies), dimnames =
   list(paste('Site', site, sep = ''), paste('Survey', survey, sep = ''),
    paste('Year', year, sep = ''), paste('Spec', spec, sep = '')))# Det. hist and p

# Create covariates (same for all species)
Xpsi1 <- matrix(runif(nsites, -2, 2), ncol = 1, dimnames = list(paste('Site',site,sep=''), NULL))  # Site covariate for psi1
Xphi <- array(runif(nsites*nyears, -2, 2), dim = c(nsites,nyears), dimnames =
 list(paste('Site',site,sep=''), paste('Year',year,sep =''))) # Yearly-site cov
Xgamma <- array(runif(nsites*nyears, -2, 2), dim = c(nsites,nyears), dimnames =
 list(paste('Site',site,sep=''), paste('Year',year,sep =''))) # Yearly-site cov
Xp <- array(runif(nsites*nsurveys*nyears,-2,2),dim=c(nsites,nsurveys,nyears), dimnames =
   list(paste('Site', site, sep = ''), paste('Survey', survey, sep = ''),
    paste('Year', year, sep = '')))  # Observation cov.

# (1) Simulate all parameter values
# (a) State process parameters
# initial occupancy for all species
mu.lpsi1 <- ifelse(mean.psi1 == '1', 500, qlogis(mean.psi1))
beta0.lpsi <- rnorm(nspecies, mu.lpsi1, sig.lpsi1)     # initial occupancy intercept
beta1.lpsi <- rnorm(nspecies, mu.beta.lpsi1, sig.beta.lpsi1) # occ. slope on Xpsi1
for(s in 1:nspecies){
  psi[,1,s] <- plogis(beta0.lpsi[s] + beta1.lpsi[s] * Xpsi1)    # psi1
}
# persistence and colonization for all species
beta0.lphi <- beta0.lgamma <- array(dim = c(nspecies, nyears-1))
mean.phi <- runif(n = nyears-1, min = min(range.mean.phi), max = max(range.mean.phi))
mean.gamma <- runif(n = nyears-1, min = min(range.mean.gamma), max = max(range.mean.gamma))
mu.lphi <- ifelse(mean.phi == '1', 500, qlogis(mean.phi))
mu.lgamma <- ifelse(mean.gamma == '1', 500, qlogis(mean.gamma))
eps.lphi <- rnorm(nspecies, 0, sig.lphi) # species effect in logit(phi) intercept
eps.lgamma <- rnorm(nspecies, 0, sig.lgamma) # spec effect in logit(gam) intercept
for(t in 1:(nyears-1)){
  beta0.lphi[,t] <- mu.lphi[t] + eps.lphi       # logit(phi) intercept
  beta0.lgamma[,t] <- mu.lgamma[t] + eps.lgamma # logit(gamma) intercept
}
beta1.lphi <- rnorm(nspecies, mu.beta.lphi, sig.beta.lphi) # slope of logit(phi) on Xphi
beta1.lgamma <- rnorm(nspecies, mu.beta.lgamma, sig.beta.lgamma) # slope of logit(gamma) on Xphi
for(s in 1:nspecies){
  for(t in 1:(nyears-1)){
    phi[,t, s] <- plogis(beta0.lphi[s, t] + beta1.lphi[s] * Xphi[,t])
    gamma[,t,s] <- plogis(beta0.lgamma[s, t] + beta1.lgamma[s] * Xgamma[,t])
  }
}

# (b) Observation process parameters
beta0.lp <- array(dim = c(nspecies, nyears))
mean.p <- runif(n = nyears, min = min(range.mean.p), max = max(range.mean.p))
mu.lp <- ifelse(mean.p == '1', 500, qlogis(mean.p))
eps.lp <- rnorm(nspecies, 0, sig.lp) # species effect in logit(p) intercept
for(t in 1:nyears){
  beta0.lp[,t] <- mu.lp[t] + eps.lp       # logit(p) intercept
}
beta1.lp <- rnorm(nspecies, mu.beta.lp, sig.beta.lp) # slope of logit(p) on Xp
beta1 <- runif(n = nyears, min = min(range.beta1.survey), max = max(range.beta1.survey))
beta2 <- runif(n = nyears, min = min(range.beta2.survey), max = max(range.beta2.survey))
sd.site <- seq(from = trend.sd.site[1], to = trend.sd.site[2], length.out = nyears)
sd.survey <- seq(from = trend.sd.survey[1], to = trend.sd.survey[2], length.out = nyears)

# Create site and survey random effects
for(i in 1:nsites){
  for(t in 1:nyears){
    eps1 <- rnorm(n = nsites, mean = 0, sd = sd.site[t])  # Site random eff.
    eps2 <- rnorm(n = nsurveys, mean = 0, sd = sd.survey[t]) # Survey random eff.
  }
}
for(s in 1:nspecies){
  for(t in 1:nyears){   # Years
    for(j in 1:nsurveys){ # Occasions interpreted as surveys
      p[,j,t,s] <- plogis(beta0.lp[s, t] + beta1.lp[s] * Xp[,j,t] +
      eps1 + eps2[j] + beta1[t] * (j - (nsurveys/2)) +
      beta2[t] * (j - (nsurveys/2))^2)
    }
  }
}

# (2) Simulate the true system dynamics (state process)
# First year
for(s in 1:nspecies){
  z[,1, s] <- rbinom(nsites, 1, psi[,1,s])   # Initial occurrence state
}

for(s in 1:nspecies){                    # Loop over species
  for(t in 2:nyears){                # Loop over years
    muZ[,t,s] <- z[,t-1,s] * phi[,t-1,s] + (1-z[,t-1,s]) * gamma[,t-1,s]
    z[,t,s] <- rbinom(nsites, 1, muZ[,t,s])
  }
}

# (3) Simulate observation process to get the observed data
for(s in 1:nspecies){                    # Loop over species
  for(t in 1:nyears){                # Loop over years
    for(j in 1:nsurveys){               # Loop over replicates
      prob <- z[,t,s] * p[,j,t,s] # zero out p for unoccupied sites
      y[,j,t,s] <- rbinom(nsites, 1, prob)
    }
  }
}

# (4) Compute annual population occupancy
for(s in 1:nspecies){                    # Loop over species
  for (t in 2:nyears){
    psi[,t,s] <- psi[,t-1,s] * phi[,t-1,s] + (1-psi[,t-1,s]) * gamma[,t-1,s]
  }
}

# Compute some derived stuff
n.occ <- apply(z, 2:3, sum)             # Number of occupied sites
psi.fs <- apply(z, 2:3, mean)           # Finite-sample occupancy proportion
mean.psi <- apply(psi, 2:3, mean)       # Average psi over sites
z.obs <- apply(y, c(1,3,4), max)        # Observed value of z matrix
n.occ.obs <- apply(z.obs, 2:3, sum)     # Observed number of occupied sites
psi.obs <- apply(z.obs, 2:3, mean)      # Observed occupancy (finite sample)

# Total number of species that occur in the sampled sites
tmp1 <- apply(z, 2:3, max)              # True presence per year and species
nyears.pres <- apply(tmp1, 2, sum)       # Number of years when species present
nspecies.pres <- sum(nyears.pres > 0)       # Number of species ever present

# Total number of species that were detected anywhere in the sampled sites
tmp2 <- apply(z.obs, 2:3, max)          # Observed presence per year and species
nyears.det <- apply(tmp2, 2, sum)        # Number of years when species detected
nspecies.det <- sum(nyears.det > 0)         # Number of species ever detected


# Print out number of occurring and detected species
if(verbose) {
  cat(paste("\n *** Number of species ever occurring:", nspecies.pres,
    "\n *** Number of species ever detected:", nspecies.det,
    "\n *** Avg. number of years of occurrence:", round(mean(nyears.pres), 3),
    "\n *** Avg. number of years with detection:", round(mean(nyears.det), 3), "\n\n"))
}
# Compute the average survey product of availability and detection
# (ignoring the other terms in the model for detection)
p.survey <- array(NA, dim = c(nsurveys, nyears))
for(t in 1:nyears){   # Years
  p.survey[,t] <- plogis(mean(beta0.lp[, t]) + beta1[t] * (survey - (nsurveys/2)) + beta2[t] * (survey - (nsurveys/2))^2)
}

# (5) Plots of stuff
if(show.plot){
  oldpar <- par(mfrow = c(3, 2), mar = c(5,5,4,3), cex.lab = 1.2)
  oldAsk <- devAskNewPage(ask = dev.interactive(orNone = TRUE))
  on.exit({par(oldpar) ; devAskNewPage(oldAsk)})

  tryPlot <- try( {
    # Get predicted covariate relationships and plot them in single graph
    pred.cov <- seq(-2, 2, length.out = 100)
    psi.pred <- phi.pred <- gamma.pred <- p.pred <- array(dim = c(length(pred.cov), nspecies))
    for(s in 1:nspecies){
      psi.pred[,s] <- plogis(beta0.lpsi[s] + beta1.lpsi[s] * pred.cov)
      phi.pred[,s] <- plogis(mean(beta0.lphi[s,]) + beta1.lphi[s] * pred.cov)
      gamma.pred[,s] <- plogis(mean(beta0.lgamma[s,]) + beta1.lgamma[s] * pred.cov)
      p.pred[,s] <- plogis(mean(beta0.lp[s,]) + beta1.lp[s] * pred.cov)
    }
    matplot(pred.cov, psi.pred, type = 'l', lty = 1, ylim = c(0,1), lwd = 2,
        main = paste('Occupancy (', nspecies, ' species, ', nsites, ' sites)', sep = ''),
        xlab = 'Covariate', ylab = 'Initial occupancy prob.', las = 1, frame = FALSE)
    matplot(pred.cov, phi.pred, type = 'l', lty = 1, ylim = c(0,1), lwd = 2,
        main = paste('Persistence (averaged over years,\n', nspecies, ' species, ', nsites, ' sites)', sep = ''), 
        xlab = 'Covariate', ylab = 'Persistence prob.', las = 1, frame = FALSE)
    matplot(pred.cov, gamma.pred, type = 'l', lty = 1, ylim = c(0,1), lwd = 2, 
        main = paste('Colonization (averaged over years,\n', nspecies, ' species, ', nsites, ' sites)', sep = ''), 
        xlab = 'Covariate', ylab = 'Colonization prob.', las = 1, frame = FALSE)
    matplot(pred.cov, p.pred, type = 'l', lty = 1, ylim = c(0,1), lwd = 2,
        main = paste('Detection (averaged over years,\n', nspecies, ' species, ', nsites, ' sites)', sep = ''), 
        xlab = 'Covariate', ylab = 'Detection prob.', las = 1, frame = FALSE)

    # Plot the average surveyal product of availability and detection
    # (ignoring the other terms in the model for detection)
    matplot(survey, p.survey, type = 'l', lty = 1, lwd = 2,
        main = 'Seasonal pattern in p over the years \n(only survey terms, same for all species)', 
        xlab = 'Survey', ylab = 'Detection probability', ylim = c(0,1))

    # Histo of detection
    hist(p, col = 'grey', breaks = 50, xlim = c(0,1), 
        main = 'Detection probability p\n (all species, sites etc.)')


    # Annual (and species-specific) variation in persistence, colonisation, and detection
    matplot(t(plogis(beta0.lphi)), type = 'l', lty = 1, lwd = 2, ylim = c(0,1), xlab = 'Year',
        ylab = 'Persistence intercept', main = 'Average persistence per year and species',
        las = 1, frame = FALSE)
    matplot(t(plogis(beta0.lgamma)), type = 'l', lty = 1, lwd = 2, ylim = c(0,1), xlab = 'Year',
        ylab = 'Colonization intercept', main = 'Average colonization per year and species',
        las = 1, frame = FALSE)
    matplot(t(plogis(beta0.lp)), type = 'l', lty = 1, lwd = 2, ylim = c(0,1), xlab = 'Year',
        ylab = 'Detection intercept', main = 'Average detection per year and species', 
        las = 1, frame = FALSE)

    # Histo of true mean occupancy probability (all species and years)
    hist(mean.psi, col = 'grey', breaks = 50, xlim = c(0,1),
        main = 'Mean occupancy probability psi1\n (all species and years)')

    # Plot realised and apparent proportion of occupied sites
    matplot(year, mean.psi, type = "l", lty = 1, xlab = "Year", ylab = "Occupancy prob.",
        xlim = c(0,nyears+1), ylim = c(0,1), lwd = 2, frame.plot = FALSE, las = 1,
        main = paste('True occupancy (', nspecies, ' species, ', nsites, ' sites)', sep = '') )
    matplot(year, psi.obs, type = "l", lty = 1, xlab = "Year", ylab = "Occupancy prob.",
        xlim = c(0,nyears+1), ylim = c(0,1), lwd = 2, frame.plot = FALSE, las = 1, 
        main = paste('Observed occupancy (', nspecies, ' species, ', nsites, ' sites)', sep = ''))
  }, silent = TRUE)
  if(inherits(tryPlot, "try-error"))
    tryPlotError(tryPlot)
}

# Return data
return(list(nspecies = nspecies, nsites = nsites, nsurveys = nsurveys, nyears = nyears, mean.psi1 = mean.psi1, sig.lpsi1 = sig.lpsi1, mu.beta.lpsi1 = mu.beta.lpsi1, sig.beta.lpsi1 = sig.beta.lpsi1, range.mean.phi = range.mean.phi, sig.lphi = sig.lphi, mu.beta.lphi = mu.beta.lphi, sig.beta.lphi = sig.beta.lphi, range.mean.gamma = range.mean.gamma, sig.lgamma = sig.lgamma, mu.beta.lgamma = mu.beta.lgamma, sig.beta.lgamma = sig.beta.lgamma, range.mean.p = range.mean.p, sig.lp = sig.lp, mu.beta.lp = mu.beta.lp, sig.beta.lp = sig.beta.lp, range.beta1.survey = range.beta1.survey, range.beta2.survey = range.beta2.survey, trend.sd.site = trend.sd.site, trend.sd.survey = trend.sd.survey, Xpsi1 = Xpsi1, Xphi = Xphi, Xgamma = Xgamma, Xp = Xp, beta0.lpsi = beta0.lpsi, beta1.lpsi = beta1.lpsi, psi = psi, mean.phi = mean.phi, mean.gamma = mean.gamma, eps.lphi = eps.lphi, eps.lgamma = eps.lgamma, beta0.lphi = beta0.lphi, beta0.lgamma = beta0.lgamma, beta1.lphi = beta1.lphi, beta1.lgamma = beta1.lgamma, phi = phi, gamma = gamma, mean.p = mean.p, eps.lp = eps.lp, beta0.lp = beta0.lp, beta1.lp = beta1.lp, beta1 = beta1, beta2 = beta2, sd.site = sd.site, sd.survey = sd.survey, eps1 = eps1, eps2 = eps2, n.occ = n.occ, psi.fs = psi.fs, mean.psi = mean.psi, z.obs = z.obs, n.occ.obs = n.occ.obs, psi.obs = psi.obs, nyears.pres = nyears.pres, nspecies.pres = nspecies.pres, nyears.det = nyears.det, nspecies.det = nspecies.det, z = z, p = p, y = y))
} # ------------------ End function definition ---------------------

