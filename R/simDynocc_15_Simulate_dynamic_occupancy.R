simDynocc <- function(nsite = 250, nrep = 3, nyear = 10, year.impact = 5,
   mean.psi1 = 0.4, beta.Xpsi1 = 0, 
   range.phi = c(0.5, 1), impact.phi = 0, beta.Xphi = 0, 
   range.gamma = c(0, 0.5), impact.gamma = 0, beta.Xgamma = 0, 
   range.p = c(0.1, 0.9), beta.Xp = 0,
   range.beta1.season = c(0, 0), range.beta2.season = c(0, 0), 
   range.sd.site = c(0, 0), range.sd.survey = c(0, 0), 
   range.sd.site.survey = c(0, 0), show.plot = TRUE) {
#
# Written by Marc Kéry, 4 Dec 2014
#
# Adapted on 18-20 October 2016
#
# Function to simulate detection/nondetection data under a general 
#     dynamic site-occ model, including:
#   * annual variation in the probabilities of patch persistence, colonization
#     and detection is specified by the bounds of a uniform distribution.
#   * one covariate is allowed to affect a parameter: a site covariate for psi1,
#     a site-by-year covariate for phi and gamma and an 
#     observational covariate for p
#   * Additional detection heterogeneity at the site-, the survey, or the 
#     site-by-survey level, with the possibility of a temporal trend in 
#     this heterogeneity over the years. E.g., an annual trend in 
#     detection heterogeneity at the site or the survey level is specified 
#     by the first and second value, which correspond to the heterogeneity in
#     the first and the last year Hence, range.sd.site = c(0, 1) will result in
#     a linear trend in the magnitude of site heterogeneity in detection 
#     from 0 in the first year to 1 in the last year.
#   * Additional detection heterogeneity that varies over the season 
#     (= occasion) according to a quadratic effect of occasion number
#     (to model phenology of an insect species for instance).
#   * Simulation of data under a BACI (before-after-control-impact) design,
#     where some event happens in a given year and *reduces* phi or gamma 
#     by a stated percentage (only reductions, no increases allowed !)

# Function arguments:
# -------------------
# ** Sample size arguments **
# nsite – Number of sites
# nrep – Number of replicate surveys within a year (= season)
# nyear – Number of years (or 'seasons')
# year.impact – Year when some impact happens (for BACI design)
#
# ** Arguments to set intercepts of regressions **
# mean.psi1 – average occupancy probability in first year
# range.p – bounds of uniform distribution from which annual p drawn 
# range.psi and range.gamma – same for survival and colonization probability
# -------------------
#
# ** Arguments to set slopes of regressions **
# beta.Xpsi1, beta.Xphi, beta.Xgamma, beta.Xp – covariate coefficients of 
#      prob. of initial occupancy, persistence, colonization and detection.
# -------------------
#
# *** Args. for detection heterogeneity among sites, surveys, and occasions 
# range.sd.site: sd of normal distribution to model logit-normal noise in p
#      at the site level in the first and the last year of the simulation
# range.sd.survey: sd of normal distribution to model logit-normal noise in p
#      ONLY at the rep = ‘survey’ level, in the first and the last year.
# range.sd.site.survey: sd of normal distribution to model logit-normal noise 
#      in p at the site/year/rep = ‘survey’ level, in the first and the 
#      last year.
# For these arguments, if the two values in the range are the
#      same, a constant value is assumed over time, while if they are different,
#      a linear trend is assumed over time.
# range.beta1.season is the range of the annual variation in the linear effect 
#     of season (i.e., of month 1-12) on the product of 
#     availability and detection linear and quadratic effect of season 
# range.beta2.season is the same for the quadratic effect of season
# -------------------
#
# ** Arguments for the BACI design **
# year.impact: year in which an impact happens, which affects phi and gamma
# impact.phi: effect in percent on annual phi (must be zero or negative, 
#     e.g., impact.phi = -20 for a 20% reduction in phi)
# impact.gamma: effect in percent on annual gamma
#

# Set up arrays needed
site <- 1:nsite                             # Sites
year <- 1:nyear                             # Years
month <- 1:nrep                            # Months (= visit)
psi <- muZ <- z <- array(dim = c(nsite, nyear), dimnames = 
list(paste('Site', site, sep = ''), paste('Year', year, sep = ''))) # Occupancy, occurrence
phi <- gamma <- array(NA, dim = c(nsite, (nyear-1)), dimnames =
   list(paste('Site', site, sep = ''), paste('Year', year[-nyear], sep = ''))) # Survival, colonisation
y <- p <- array(NA, dim = c(nsite, nrep, nyear), dimnames =
   list(paste('Site', site, sep = ''), paste('Visit', month, sep = ''),
    paste('Year', year, sep = '')))# Det. hist and p

# Create covariates
Xpsi1 <- runif(nsite, -2, 2)                     # Site covariate for psi1
Xphi <- array(runif(nsite*nyear, -2, 2), dim = c(nsite,nyear))
                                             # Yearly-site covariate
Xgamma <- array(runif(nsite*nyear, -2, 2), dim = c(nsite,nyear))
                                             # Yearly-site covariate
Xp <- array(runif(nsite*nrep*nyear,-2,2),dim=c(nsite,nrep,nyear))  # Observ. cov.

# Create impact covariate and effect on phi and gamma
impact <- rep(0, (nyear-1))
impact[year.impact:(nyear-1)] <- 1

# (1) Simulate all parameter values
# (a) State process parameters
psi[,1] <- plogis(qlogis(mean.psi1) + beta.Xpsi1 * Xpsi1) # psi1
mean.phi <- runif(n = nyear-1, min = range.phi[1], max = range.phi[2])
mean.gamma <- runif(n = nyear-1, min = range.gamma[1], max = range.gamma[2])

# effect of impact on phi: additive effect on persistence
phi.effect <- (mean.phi * (impact.phi/100) * impact)
# effect of impact on gamma: additive effect on colonisation
gamma.effect <- (mean.gamma * (impact.gamma/100) * impact)

# Assemble effects of year, impact and covariates on phi and gamma
for(t in 1:(nyear-1)){
   phi[,t] <- plogis(qlogis(mean.phi[t] + phi.effect[t]) + beta.Xphi * Xphi[,t])
   gamma[,t] <- plogis(qlogis(mean.gamma[t] + gamma.effect[t]) + beta.Xgamma * Xgamma[,t])
}


# (b) Observation process parameters
mean.p <- runif(n = nyear, min = range.p[1], max = range.p[2])
beta1 <- runif(n = nyear, min = range.beta1.season[1], max = range.beta1.season[2])
beta2 <- runif(n = nyear, min = min(range.beta2.season), max = max(range.beta2.season))
# Next two allow incorporation of trend over time
sd.site <- seq(from = range.sd.site[1], to = range.sd.site[2], length.out = nyear)
sd.survey <- seq(from = range.sd.survey[1], to = range.sd.survey[2], length.out = nyear)
sd.site.survey <- seq(from = range.sd.site.survey[1], to = range.sd.site.survey[2], length.out = nyear)

# Define and fill the array of site/survey random effects
eps3 <- array(dim = c(nsite, nrep, nyear))
for(t in 1:nyear){
  eps3[,,t] <- matrix(rnorm(n = nsite*nrep, sd = sd.site.survey[t]), ncol = nrep)
}

for(i in 1:nsite){     # Sites
  for(t in 1:nyear){   # Years
    eps1 <- rnorm(n = nsite, sd = sd.site[t])   # Zero-mean site random eff.
    eps2 <- rnorm(n = nrep, sd = sd.survey[t]) # Zero-mean survey random eff.
      # ZM site.survey ranef.
    for(j in 1:nrep){ # Months
      p[i,j,t] <- plogis(qlogis(mean.p[t]) + beta.Xp*Xp[i,j,t] + 
        eps1[i] + eps2[j] + eps3[i,j,t] + 
        beta1[t] * (j - (nrep/2)) + beta2[t] * (j - (nrep/2))^2)
    }
  }
}


# (2) Simulate the true system dynamics (state process)
# First year
z[,1] <- rbinom(nsite, 1, psi[,1])   # Initial occurrence state
# Years 2:nyear
for(i in 1:nsite){                   # Loop over sites
   for(t in 2:nyear){                # Loop over years
      muZ[i,t] <- z[i, t-1]*phi[i,t-1] + (1-z[i, t-1])*gamma[i,t-1]
      z[i,t] <- rbinom(1, 1, muZ[i,t])
   }
}

# (3) Simulate observation process to get the observed data
for(i in 1:nsite){                     # Loop over sites
   for(t in 1:nyear){                  # Loop over years
      for(j in 1:nrep){               # Loop over replicates
         prob <- z[i,t] * p[i,j,t] # zero out p for unoccupied sites
         y[i,j,t] <- rbinom(1, 1, prob)
      }
   }
}

# (4) Compute annual population occupancy
for(i in 1:nsite){
   for (t in 2:nyear){
      psi[i,t] <- psi[i,t-1]*phi[i,t-1] + (1-psi[i,t-1])*gamma[i,t-1]
   }
}
n.occ <- apply(z, 2, sum)             # Number of occupied sites
psi.fs <- apply(z, 2, mean)           # Finite-sample occupancy proportion
mean.psi <- apply(psi, 2, mean)                   # Average psi over sites
psi.app <- apply(apply(y, c(1,3), max), 2, mean)  # Apparent occupancy (finite sample)

# Compute average seasonal product of availability and detection
# (ignoring the other terms in the model for detection)
p.season <- array(NA, dim = c(nrep, nyear))
for(t in 1:nyear){   # Years
  p.season[,t] <- plogis(qlogis(mean.p[t]) + beta1[t] * (month - (nrep/2)) + beta2[t] * (month - (nrep/2))^2)
}


# (5) Plots of stuff
if(show.plot){
  op <- par(mfrow = c(2, 2)) ; on.exit(par(op))

  # Get predicted covariate relationships and plot them in single graph
  pred.cov <- seq(-2, 2, length.out = 100)
  psi.pred <- plogis(qlogis(mean.psi1) + beta.Xpsi1 * pred.cov)
  phi.pred <- plogis(mean(qlogis(mean.phi)) + beta.Xphi * pred.cov)
  gamma.pred <- plogis(mean(qlogis(mean.gamma)) + beta.Xgamma * pred.cov)
  p.pred <- plogis(mean(qlogis(mean.p)) + beta.Xp * pred.cov)

  plot(pred.cov, psi.pred, type = 'l', col = 'green', ylim = c(0,1), lwd = 2, main = 'Covariate relationships\n (green - psi1, blue - phi, black - gamma, red - p)', xlab = 'Covariate value', ylab = 'Predicted Prob.')
  lines(pred.cov, phi.pred, type = 'l', col = 'blue', lwd = 2)
  lines(pred.cov, gamma.pred, type = 'l', col = 'black', lwd = 2)
  lines(pred.cov, p.pred, type = 'l', col = 'red', lwd = 2)

  # Seasonal pattern of detection (= product of availability and detection)
  # (ignoring the other terms in the model for detection)
  matplot(month, p.season, type = 'l', lty = 1, lwd = 2, main = 'Seasonal pattern in p over the years \n(only seasonal terms)', xlab = 'Month', ylab = 'Detection probability', ylim = c(0,1))

  # Histo of detection
  hist(p, col = 'grey', breaks = 100, xlim = c(0,1), main = 'Detection probability p')

  # Plot realised and apparent proportion of occupied sites
  plot(year, apply(z, 2, mean), type = "l", xlab = "Year", ylab = "Occupancy / Detection prob.", col = "red", xlim = c(0,nyear+1), ylim = c(0,1), lwd = 2, lty = 1, frame.plot = FALSE, las = 1, main = 'True and observed states and annual p')
  lines(year, mean.p , type = "l", col = "red", lwd = 2, lty = 2)
  lines(year, psi.app, type = "l", col = "black", lwd = 2)
  text(0.3*nyear, 0.9, labels = "red solid - true (finite sample) occupancy\n red dashed - detection\n black - observed occupancy", cex = 0.8)
}

# Return data
return(list(nsite=nsite, nrep=nrep, nyear=nyear, year.impact = year.impact, impact = impact, mean.psi1=mean.psi1, beta.Xpsi1=beta.Xpsi1,
range.phi=range.phi, impact.phi = impact.phi, beta.Xphi=beta.Xphi, phi.effect = phi.effect, range.gamma=range.gamma, impact.gamma = impact.gamma, beta.Xgamma=beta.Xgamma, gamma.effect = gamma.effect, range.p=range.p, beta.Xp=beta.Xp, range.sd.site=range.sd.site, range.sd.survey=range.sd.survey, 
range.beta1.season = range.beta1.season, range.beta2.season = range.beta2.season, beta1 = beta1, beta2 = beta2, p.season = p.season,
sd.site=sd.site, sd.survey=sd.survey, mean.phi=mean.phi, mean.gamma=mean.gamma, mean.p=mean.p, psi=psi, mean.psi=mean.psi, n.occ = n.occ, psi.fs = psi.fs, psi.app=psi.app, z=z, phi=phi, gamma=gamma, p=p, y = y, Xpsi1 = Xpsi1, Xphi = Xphi, Xgamma = Xgamma, Xp = Xp, eps3 = eps3))
} 


