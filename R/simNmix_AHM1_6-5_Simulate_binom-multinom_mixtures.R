# Functions for the book Applied Hierarchical Modeling in Ecology (AHM)
# Marc Kery & Andy Royle, Academic Press, 2016.

# simNmix - AHM1 section 6.5 p241

# Function to simulate data for binomial and multinomial mixture models under wide range of conditions (introduced in AHM1 Section 6.5)
simNmix <- function(nsites = 267, nvisits = 3, mean.theta = 1, mean.lam = 2, mean.p = 0.6, area = FALSE, beta1.theta = 0, beta2.theta = 0, beta3.theta = 0, beta2.lam = 0, beta3.lam = 0, beta4.lam = 0, beta3.p = 0, beta5.p = 0, beta6.p = 0, beta.p.survey = 0, beta.p.N = 0, sigma.lam = 0, dispersion = 10, sigma.p.site = 0, sigma.p.visit = 0, sigma.p.survey = 0, sigma.p.ind = 0, Neg.Bin = FALSE, open.N = FALSE, show.plots = TRUE, verbose = TRUE) {
#
# This very general function generates single-season count data
# under variants of the binomial N-mixture model of Royle (2004) and of the
# multinomial N-mixture model of Royle et al (2007).
#
# Data are simulated at the level of each individual and individual-specific
#   detection heterogeneity can be included. As a side-effect, individual-
#   specific detection histories are generated and hence, data are also
#   be simulated under the corresponding multinomial N-mixture model.
#
# Broadly, the function can generate data under this most general model:
#
# 'Suitability' (zero-inflation) ~ cov1 + cov2 + cov3
#
# Abundance ~ offset + cov2 + cov3 + cov4 + overdispersion
#
# Detection ~ cov3 + cov5 + cov6 + survey.covariate +
#   log(N+1) + eps.site + eps.visit + eps.survey + eps.individual
#
# Overdispersion in abundance is modelled either as a Poisson-log-normal with
#   a normal random site effect in lambda or with a Negative binomial with
#   mean lambda and a 'size', or dispersion, parameter.
#   Variable site areas can be specified to affect abundance as in an offset.
# Abundance can be zero-inflated (this is the 'suitability' model). Note that
#   the zero-inflation parameter is called theta here (in unmarked it is called
#   psi). mean.phi is the probability that a site is suitable (i.e., 1 minus
#   the expected proportion of sites with structural zero abundance.
# Site covariate 2 can affect both suitability and abundance, while covariate 3
#   may affect all three levels. Hence, the function permits to simulate the
#   case where a single site covariate affects different levels in the process
#   (e.g., abundance and detection) in opposing directions (as for instance
#   in Kery, Auk, 2008)
# Density-dependent detection can be modelled as a logistic-linear effect
#   of local abundance (centered and log(x+1) transformed)
# Overdispersion in detection is modelled via normal random effects (the eps
#   terms above) specific to sites, visits, surveys or individuals.
# Effects of covariates and random-effects factors are modelled
#   as additive on the link scale (log for abundance and
#   logit for suitability and detection).
#
# Data may be generated under one specific open-population model when
#    argument 'open.N' is set to TRUE.
#
# Written by Marc Kery, 2014-2015
#
# Function arguments
# nsites: number of sites
# nvisits: number of visits per site
# mean.theta: proportion of sites that can have non-zero abundance in principle:
#   suitability model for zero-inflation
# mean.lam: Expected abundance at the average value of all
#   abundance covariates (and ignoring random site effects): abundance model
# mean.p: Expected detection at the average value of all
#   detection covariates (and ignoring all random effects): detection model
# area: determines area of sites (A), defaults to A=1 (i.e., all identical),
#   but you can supply a vector of site areas of length nsites instead.
# beta1.theta: coefficient of site covariate 1 in suitability model
# beta2.theta: coefficient of site covariate 2 in suitability model
# beta3.theta: coefficient of site covariate 3 in suitability model
# beta2.lam: coefficient of site covariate 2 in abundance model
# beta3.lam: coefficient of site covariate 3 in abundance model
# beta4.lam: coefficient of site covariate 4 in abundance model
# beta3.p: coefficient of site covariate 3 in detection model
# beta5.p: coefficient of site covariate 5 in detection model
# beta6.p: coefficient of site covariate 6 in detection model
# beta.p.survey: coefficient of survey ('observational') covariate on p
# beta.p.N: coefficient of centered local population size (log(N+1)) in
#    detection model (i.e., coef. for density-dependent detection prob.)
# sigma.lam: "Overdispersion SD" in linear predictor of abundance
# dispersion: 'size' or extra-Poisson dispersion of Negative binomial
# sigma.p.site: "Overdispersion SD" in linear predictor of
#   detection coming from random site effects
# sigma.p.visit: "Overdispersion SD" in linear predictor of
#   detection coming from random visit effects
# sigma.p.survey: "Overdispersion SD" in linear predictor of
#   detection coming from random site-by-survey effects
# sigma.p.ind: "Overdispersion SD" in linear predictor of
#   detection coming from random site-by-individual effects
# Neg.Bin: if FALSE, any overdispersion in abundance is modelled by
#   a Poisson log-normal; if TRUE, abundance overdispersion is modelled
#   by adoption of a Negative binomial distribution for latent N
# Open.N: if TRUE, data are simulated under one specific form of an open
#   population, where N in the first occasion is drawn from the specified
#   mixture distribution and for all further occasions j, we have
#   N_ij ~ Poisson(N_i(j-1)). With open.N = TRUE, we must have
#   sigma.p.ind = 0, show.plots = FALSE and nvisits >1.
# show.plots: if TRUE, plots of the data will be displayed; set to FALSE
#      if you are running simulations.

if(FALSE) x <- NULL # Fix issues with 'curve'
logit <- plogis # allows 'logit' to appear in axis label instead of 'plogis'

# Checks and fixes for input data -----------------------------
nsites <- round(nsites[1])
nvisits <- round(nvisits[1])
stopifnotProbability(mean.theta)
stopifNegative(mean.lam, allowZero=FALSE)
stopifnotProbability(mean.p)
stopifNegative(sigma.lam)
stopifNegative(dispersion, allowZero=FALSE)
stopifNegative(sigma.p.site)
stopifNegative(sigma.p.visit)
stopifNegative(sigma.p.survey)
stopifNegative(sigma.p.ind)
# --------------------------------------------

# Create indices
nreps <- rep(nvisits, nsites)                   # No. visits (reps) per site
site <- 1:nsites                              # Site index at site level
site.per.unit <- rep(1:nsites, each = nvisits) # Site index at rep level

if(verbose) {
  cat("***** New simulation *****\n\n")
  cat("No. sites visited:      ", nsites, "\n")
  cat("No. rep. visits:          ", nvisits, "\n")
  cat("Total no. visits:       ", sum(nreps), "\n\n")
}
# Generate covariates with standardised values between -2 and 2
# Site covariates 1-6
site.cov <- matrix(runif(n = nsites*6, -2, 2), ncol = 6)
colnames(site.cov) <- c("cov1", "cov2", "cov3", "cov4", "cov5", "cov6")
# Survey covariate
survey.cov <- matrix(runif(n = nsites*nvisits, -2, 2), ncol = nvisits)

# get site-specific values for offset
if(area[1] == FALSE) A <- rep(1, nsites) # means no offset
if(area[1] != FALSE) A <- area          # use supplied vector as area for offset

# Simulate ecological process:
# (1) 'suitability' (leading to zero-inflation)
# Zero-inflation: create "suitability" indicator z
# Linear predictor of suitability model
alpha.theta <- ifelse(mean.theta == 1, 25, qlogis(mean.theta))   # Avoid Inf.
theta <- plogis(alpha.theta + beta1.theta * site.cov[,1] + beta2.theta * site.cov[,2] +
    beta3.theta * site.cov[,3])
s <- rbinom(n = nsites, 1, theta)  # Suitability indicator

# (2) Abundance process (for sites suitable in principle)
# Linear predictor of abundance model excluding random effects
# this is directly the lin.pred. for the Neg.Bin abundance model
log.lam.partial <- log(A) + log(mean.lam) + beta2.lam * site.cov[,2] +
    beta3.lam * site.cov[,3] + beta4.lam * site.cov[,4]

# Draw abundance under the (zero-inflated) Poisson distribution
# (For baseline comparison of the abundance distributions in histo below)
N.P <- rpois(n = nsites, lambda = s * exp(log.lam.partial))

# Draw abundance under the (zero-inflated) negative binomial distribution
N.NB <- rnbinom(n = nsites, mu = s * exp(log.lam.partial), size = dispersion)

# Draw abundance under the (zero-inflated) Poisson log-normal distribution
# Random site effects in lambda, zero out if Neg.Bin == TRUE
eta.lam <- rnorm(n = nsites, sd = sigma.lam * (1 - Neg.Bin))
# Linear predictor of PLN abundance model including random effects
log.lam <- log.lam.partial + eta.lam
# Draw realised values of abundance at each site
N.PLN <- rpois(n = nsites, lambda = s * exp(log.lam))

if(Neg.Bin == TRUE){   # Negative-binomial N's fed into variable N ....
  N <- N.NB
} else { # ... or else those from PLN mixture
  N <- N.PLN
}
Ntotal <- sum(N)     # Add up N over all M sites

# Ecological process when population open (open.N == TRUE)
N.open <- matrix(NA, nrow = nsites, ncol = nvisits)
if(open.N){
  N.open[,1] <- N
  for(j in 2:nvisits){
     N.open[,j] <- rpois(nsites, N.open[,j-1])
  }
  #cor(N.open)
  #matplot(1:nvisits, t(N.open), type = 'l')
}

# Visualization of suitability and abundance code now moved to line 281

# Simulate observation process conditional on true state N
# Create structures to be filled
nslice <- max(N)+1              # Max. number of slices in DH 3D array
inds <- DH <- p <- logit.p.partial <- logit.p <- array(NA, dim = c(nsites, nvisits, nslice))
C <- eta.p.survey <- array(NA, dim = c(nsites, nvisits))

# Determine occupied sites and table of 'existing' individuals (inds)
occ.sites <- which(N>0)
for(i in occ.sites){
   inds[i,,1:N[i]] <- 1
}

# Draw random site effects in p
eta.p.site <- rnorm(n=nsites, mean = 0, sd = sigma.p.site)
# Draw random visit effects in p
eta.p.visit <- rnorm(n=nvisits, mean = 0, sd = sigma.p.visit)
# Draw random survey (= site-by-survey) effects in p
eta.p.survey <- matrix(rnorm(n = nsites*nvisits, sd = sigma.p.survey), nrow = nsites, ncol = nvisits)
# Draw random individual (= site-by-ind) effects in p (NOT site-ind-visit !)
eta.p.ind <- array(rnorm(n = nsites* max(N), mean = 0, sd = sigma.p.ind), dim = c(nsites, nslice))
#eta.p.ind[N == 0,] <- NA            # NA out effects for non-existing inds.
for(i in 1:nsites){
   eta.p.ind[i,(N[i]+1):nslice] <- NA
}

# For default closed population (open.N == FALSE)
if(open.N == FALSE){
# Sample individuals to get individual detection histories (DH) for each site (note that DH[i,,] == 'NA' when N[i] = 0)
for(i in occ.sites){      # Loop over occupied sites (with N>0)
   for(j in 1:nvisits){    # Loop over visits
      for(n in 1:nslice){ # Loop over individuals
      # Linear predictor of detection model excl. random effects
      logit.p.partial[i,j,n] <- qlogis(mean.p) + beta3.p * site.cov[i,3] +
         beta5.p * site.cov[i,5] + beta6.p * site.cov[i,6] +
         beta.p.survey * survey.cov[i,j] +
         beta.p.N * (log(N[i]+1)-mean(log(N[i]+1)))
      # Linear predictor of detection model including random effects
      logit.p[i,j,n] <- logit.p.partial[i,j,n] + eta.p.site[i] +
         eta.p.visit[j] + eta.p.survey[i,j] + eta.p.ind[i,n]
      # Apply inverse link function
      p[i,j,n] <- plogis(logit.p[i,j,n])
      # Get individual detection histories: NA out non-existing inds.
      # (i.e., at sites where N=0)
      # prob <- inds[i,j,n] * p[i,j,n]     # NA out non-existing individuals
      # DH[i,j,n] <- rbinom(n=1, size = 1, prob = prob) # a warning every time prob is NA
      if(!is.na(inds[i,j,n]))
        DH[i,j,n] <- rbinom(n=1, size = 1, prob = p[i, j, n]) ## MM 2017-03-10
      # else ... the value stays as NA
    }
  }
}
# DH <- DH[,,-nslice]    # Get rid of unused last slice
# p <- p[,,-nslice]      # Get rid of unused last slice
DH <- DH[,,-nslice, drop=FALSE]    # Get rid of unused last slice ## MM 2017-03-10
p <- p[,,-nslice, drop=FALSE]      # Get rid of unused last slice

# Get counts C by tallying up detection histories (DH)
# Also get the sum over sites of max counts
# Account for possible single-visit data
if(length(dim(DH)) == 3){      # for typical multi-visit design
   C <- apply(DH, c(1,2), sum, na.rm = TRUE)
   summax <- sum(apply(C, 1, max))
}
if(length(dim(DH)) == 2){      # to account single-visit design
   C <- apply(DH, 1, sum, na.rm = TRUE)
   summax <- sum(C)
}
pp <- N.open <- NA     # Not available for open.N == FALSE
} # end open.N == FALSE


# For open population (open.N == TRUE)
if(open.N == TRUE){
pp <- matrix(NA, nrow = nsites, ncol = nvisits) # Define detection prob
for(j in 1:nvisits){    # Loop over visits
   # Full linear predictor of detection model
   pp[,j] <- plogis(qlogis(mean.p) + beta3.p * site.cov[i,3] +
         beta5.p * site.cov[i,5] + beta6.p * site.cov[i,6] +
         beta.p.survey * survey.cov[i,j] +
         beta.p.N * (log(N[i]+1)-mean(log(N[i]+1))) +
         eta.p.site + eta.p.visit[j] + eta.p.survey[,j])
   C[,j] <- rbinom(nsites, N.open[,j], pp[,j])
}
summax <- sum(apply(C, 1, max))
# Fill some things used in the function output
p <- pp  ;  DH <- NA
}

if(show.plots){
  # Restore graphical settings on exit ---------------------------
  oldpar <- par("mfrow", "cex", "cex.main")
  oldAsk <- devAskNewPage(ask = dev.interactive(orNone=TRUE))
  on.exit({par(oldpar); devAskNewPage(oldAsk)})
  # --------------------------------------------------------------

  tryPlot <- try( {
    # Visualization of suitability and abundance
    # ''''''''''''''''''''''''''''''''''''''''''
    # Page 1: Plots features of the suitability part of the system
    par(mfrow = c(2, 2), cex.main = 1)
    barplot(table(s), main = "Number unsuitable and suitable sites", col = "grey")
    plot(site.cov[,1], s, ylim = c(0,1), main = "'Suitability' & site covariate 1")
    curve(logit(alpha.theta + beta1.theta * x), -2, 2, col = "red", add = TRUE, lwd = 3)
    plot(site.cov[,2], s, ylim = c(0,1), main = "'Suitability' & site covariate 2")
    curve(logit(alpha.theta + beta2.theta * x), -2, 2, col = "red", add = TRUE, lwd = 3)
    plot(site.cov[,3], s, ylim = c(0,1), main = "'Suitability' & site covariate 3")
    curve(logit(alpha.theta + beta3.theta * x), -2, 2, col = "red", add = TRUE, lwd = 3)

    # Page 2: Plots features of the abundance part of the system
    par(mfrow = c(3, 3), cex.main = 1)
    ylim = c(min(exp(log.lam.partial))-1, max(N))
    curve(exp(log(mean.lam) + beta2.lam * x), -2, 2, xlab = "Site covariate 2", main = "Site covariate 2 & lambda", ylab = "partial lambda", col = "red", lwd = 3)
    curve(exp(log(mean.lam) + beta3.lam * x), -2, 2, xlab = "Site covariate 3", main = "Site covariate 3 & lambda", ylab = "partial lambda", col = "red", lwd = 3)
    curve(exp(log(mean.lam) + beta4.lam * x), -2, 2, xlab = "Site covariate 4", main = "Site covariate 4 & lambda", ylab = "partial lambda", col = "red", lwd = 3)
    plot(site.cov[,2], exp(log.lam.partial), col = "red", xlab = "Site covariate 2", ylab = "lambda", main = "Marginal lambda \n(excl. site random effects)", ylim = ylim)
    plot(site.cov[,3], exp(log.lam.partial), col = "red", xlab = "Site covariate 3", ylab = "lambda", main = "Marginal lambda \n(excl. site random effects)", ylim = ylim)
    plot(site.cov[,4], exp(log.lam.partial), col = "red", xlab = "Site covariate 4", ylab = "lambda", main = "Marginal lambda \n(excl. site random effects)", ylim = ylim)
    plot(site.cov[,2], exp(log.lam), col = "red", xlab = "Site covariate 2", ylab = "lambda", main = "Marginal lambda \n(incl. site random effects)", ylim = ylim)
    plot(site.cov[,3], exp(log.lam), col = "red", xlab = "Site covariate 3", ylab = "lambda", main = "Marginal lambda \n(incl. site random effects)", ylim = ylim)
    plot(site.cov[,4], exp(log.lam), col = "red", xlab = "Site covariate 4", ylab = "lambda", main = "Marginal lambda \n(incl. site random effects)", ylim = ylim)

    # Page 3: Realized adundances
    par(mfrow = c(1, 3), cex = 1)
    plot(site.cov[,2], N, col = "red", xlab = "Site covariate 2", ylab = "N", main = "Realized abundance (N)", ylim = ylim)
    plot(site.cov[,3], N, col = "red", xlab = "Site covariate 3", ylab = "N", main = "Realized abundance (N)", ylim = ylim)
    plot(site.cov[,4], N, col = "red", xlab = "Site covariate 4", ylab = "N", main = "Realized abundance (N)", ylim = ylim)

    # Page 4: Random site effects if !Neg.Bin, histogram for N for both
    if(Neg.Bin == TRUE){
      xlim <- c(min(c(N.P, N.NB)), max(c(N.P, N.NB)))
      par(mfrow = c(1, 1), cex.main = 1)
      histCount(N.P, N.NB, xlab = "Abundance N",
      main = paste("N under (zero-infl.) Neg.bin (red)", "and (zero-infl.) Poisson (blue) mixtures", sep="\n"))
    } else {
      xlim <- c(min(c(N.P, N.PLN)), max(c(N.P, N.PLN)))
      par(mfrow = c(1, 2), cex.main = 1)
      hist(eta.lam, col = "grey", main = "Random site effects in abundance")
      histCount(N.P, N.PLN, xlab = "Abundance N",
          main = paste(c("N under (zero-infl.) Poisson log-normal (red)",
          "compared with baseline (zero-infl.) Poisson mixture (blue)"), sep="\n"))
    }

    # Plots and summaries of observation process
    # ''''''''''''''''''''''''''''''''''''''''''
    # Page 5: Effects on p
    par(mfrow = c(3,2), cex.main = 1)
    curve(logit(qlogis(mean.p) + beta3.p * x), -2, 2, xlab = "Site covariate 3", main = "Site covariate 3 & detection", ylab = "p", col = "red", lwd = 3)
    curve(logit(qlogis(mean.p) + beta5.p * x), -2, 2, xlab = "Site covariate 5", main = "Site covariate 5 & detection", ylab = "p", col = "red", lwd = 3)
    curve(logit(qlogis(mean.p) + beta6.p * x), -2, 2, xlab = "Site covariate 6", main = "Site covariate 6 & detection", ylab = "p", col = "red", lwd = 3)
    curve(logit(qlogis(mean.p) + beta.p.survey * x), -2, 2, xlab = "Survey covariate", main = "Survey covariate & detection", ylab = "p", col = "red", lwd = 3)
    curve(logit(qlogis(mean.p) + beta.p.N * x), log(0+1), log(max(N)+1), xlab = "Effect of log(N+1) in logit(p)", ylab = "p", col = "red", lwd = 3)

    # Page 6: Random effects in p
    par(mfrow = c(2,2), cex.main = 1)
    hist(eta.p.site, col = "grey", main = "Random site eff. in p", breaks = 50)
    hist(eta.p.visit, col = "grey", main = "Random visit eff. in p", breaks = 50)
    hist(eta.p.survey, col = "grey", main = "Random site-survey eff. in p", breaks = 50)
    hist(eta.p.ind, col = "grey", main = "Random ind. eff. in p", breaks = 50)

    # Page 7: partial p and p (site covars)
    par(mfrow = c(3,2), cex.main = 1)
    matplot(site.cov[,3], apply(plogis(logit.p.partial), c(1,2), mean, na.rm = TRUE), col = "red", xlab = "Site covariate 3", ylab = "Partial p", main = "Partial expected detection \n(no random effects)", ylim = c(0,1), pch = 1)
    matplot(site.cov[,3], apply(p, c(1,2), mean, na.rm = TRUE), col = "red", xlab = "Site covariate 3", ylab = "p", main = "Detection probability (with random effects)", ylim = c(0,1), pch = 1)
    matplot(site.cov[,5], apply(plogis(logit.p.partial), c(1,2), mean, na.rm = TRUE), col = "red", xlab = "Site covariate 5", ylab = "Partial p", main = "Partial expected detection \n(no random effects)", ylim = c(0,1), pch = 1)
    matplot(site.cov[,5], apply(p, c(1,2), mean, na.rm = TRUE), col = "red", xlab = "Site covariate 5", ylab = "p", main = "Detection probability (with random effects)", ylim = c(0,1), pch = 1)
    matplot(site.cov[,6], apply(plogis(logit.p.partial), c(1,2), mean, na.rm = TRUE), col = "red", xlab = "Site covariate 6", ylab = "Partial p", main = "Partial expected detection \n(no random effects)", ylim = c(0,1), pch = 1)
    matplot(site.cov[,6], apply(p, c(1,2), mean, na.rm = TRUE), col = "red", xlab = "Site covariate 6", ylab = "p", main = "Detection probability (with random effects)", ylim = c(0,1), pch = 1)

    # Page 8: partial p and p (survey covars), p and realised p
    par(mfrow = c(2,2), cex.main = 1)
    matplot(survey.cov, apply(plogis(logit.p.partial), c(1,2), mean, na.rm = TRUE), col = "red", xlab = "Survey covariate", ylab = "Partial p", main = "Partial expected detection \n(no random effects)", ylim = c(0,1), pch = 1)
    matplot(survey.cov, apply(p, c(1,2), mean, na.rm = TRUE), col = "red", xlab = "Survey covariate", ylab = "p", main = "Detection probability (with random effects)", ylim = c(0,1), pch = 1)
    matplot(N, apply(p, c(1,2), mean, na.rm = TRUE), col = "red", xlab = "Local abundance (N)", ylab = "p", main = "Detection probability", ylim = c(0,1), pch = 1)
    hist(p, col = "grey", main = "Realized detection probability \n(blue=mean)", breaks = 50)
    abline(v = mean(p, na.rm = TRUE), col = "blue", lwd = 2)

    # Page 9: Observed counts
    par(mfrow = c(3,3), cex.main = 1)
    # hist(C, col = "grey", main = "Observed counts", breaks = 50)
    histCount(C, NULL, color = "grey", main = "Observed counts", xlab = "C")
    matplot(site.cov[,1], C, xlab = "Site covariate 1", ylab = "Counts", main = "Obs. counts vs. site covariate 1")
    matplot(site.cov[,2], C, xlab = "Site covariate 2", ylab = "Counts", main = "Obs. counts vs. site covariate 2")
    matplot(site.cov[,3], C, xlab = "Site covariate 3", ylab = "Counts", main = "Obs. counts vs. site covariate 3")
    matplot(site.cov[,4], C, xlab = "Site covariate 4", ylab = "Counts", main = "Obs. counts vs. site covariate 4")
    matplot(site.cov[,5], C, xlab = "Site covariate 5", ylab = "Counts", main = "Obs. counts vs. site covariate 5")
    matplot(site.cov[,6], C, xlab = "Site covariate 6", ylab = "Counts", main = "Obs. counts vs. site covariate 6")
    matplot(survey.cov, C, xlab = "Survey covariate", ylab = "Counts", main = "Obs. counts vs. survey covariate")
    plot(rep(N, nvisits), C, xlab = "True state (abundance N)", ylab = "Obs.state (counts C)", main = "Obs. counts vs. true abundance", xlim = c(min(N,C), max(N,C)), ylim = c(min(N,C), max(N,C)))
    abline(0,1)
  }, silent = TRUE )
  if(inherits(tryPlot, "try-error"))
    tryPlotError(tryPlot)
}

# Compute naive 'overdispersion coefficients' at level latent N and observed C
odcN <- round(var(N)/mean(N),2)     # Overdispersion coefficient
if(open.N){
  odcN <- round(var(c(N.open))/mean(N.open),2)     # Overdispersion coefficient
}
odcC <- round(var(c(C))/mean(C),2)     # Overdispersion coefficient
if(verbose) {
  cat("\nNaive overdispersion measure (var/mean) for true abundance (N):", odcN,"\n")
  cat("Naive overdispersion measure (var/mean) for observed counts (C):", odcC,"\n")
}

# Output
# *** Key output elements are ***
# DH: detection history for each of N individuals detected at the nsites sites
# C: summary of DH: number of individuals detected for each site and visit
#
return(list(nsites = nsites, nvisits = nvisits, nobs = sum(nreps), Neg.Bin = Neg.Bin, open.N = open.N, area = area, mean.theta = mean.theta, mean.lam = mean.lam, mean.p = mean.p, beta1.theta = beta1.theta, beta2.theta = beta2.theta, beta3.theta = beta3.theta, beta2.lam = beta2.lam, beta3.lam = beta3.lam, beta4.lam = beta4.lam, beta3.p = beta3.p, beta5.p = beta5.p, beta6.p = beta6.p, beta.p.survey = beta.p.survey, beta.p.N = beta.p.N, sigma.lam = sigma.lam, dispersion = dispersion, sigma.p.site = sigma.p.site, sigma.p.visit = sigma.p.visit, sigma.p.survey = sigma.p.survey, sigma.p.ind = sigma.p.ind, site.cov = site.cov, survey.cov = survey.cov, log.lam = log.lam, s = s, N = N, p = p, DH = DH, N.open = N.open, C = C, eta.lam = eta.lam, eta.p.site = eta.p.site, eta.p.visit = eta.p.visit, eta.p.survey = eta.p.survey, eta.p.ind = eta.p.ind, odcN = odcN, odcC = odcC, Ntotal = Ntotal, summax = summax))
}

