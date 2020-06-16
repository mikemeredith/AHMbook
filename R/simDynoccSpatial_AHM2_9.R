
# 1. Define an R function to generate dynamic presence/absence systems with 'space'

# Code to define a function for simulating data.

simDynoccSpatial <- function(side = 50, nyears = 10, nsurveys = 3,
      mean.psi1 = 0.4, beta.Xpsi1 = 0,
      range.phi = c(0.8, 0.8), beta.Xphi = 0,
      range.gamma = c(0.1, 0.1), beta.Xgamma = 0,
      range.p = c(0.4, 0.4), beta.Xp = 0,
      theta.XAC = 5000, beta.XAC = c(0, 0, 0, 0), beta.Xautolog = c(0, 0),
      trend.sd.site = c(0, 0), trend.sd.survey = c(0, 0),
      seed.XAC = NA, seed = NULL, show.plots= TRUE, ask.plot = TRUE, verbose=TRUE) {
  #
  #      Written by Marc Kéry, 2014-2018
  #
  # Function to simulate detection/nondetection data in a square area
  #   under a very general dynamic site-occ model, including the
  #   following effects:
  # (1) annual variation in the probabilities of patch persistence,
  #   colonization and detection can be specified by the bounds of a
  #   uniform distribution.
  # (2) one site-, site/year-, and site/year/rep-specific covariate
  #   is allowed to affect the probabilities of occupancy
  #   (beta.Xpsi1 for site-covariate), colonisation/persistence
  #   (beta.Xgamma, beta.Xphi, for yearly site-covariate), and
  #   detection (beta.Xp for observational covariate), respectively.
  # (3) a single, spatially structured covariate for habitat suitability
  #   may affect all parameters via coefficient beta.XAC (for a
  #   biologically reasonable way, choose coefficients with the same sign
  #   for all 4 (mediated by underlying density).
  #   That spatial covariate is simulated as a Gaussian random field
  #   with negative exponential correlation function with
  #   'range parameter' theta.XAC
  # (4) autologistic effects (beta.Xautolog) in persistence and colonization
  #   probability can be chosen, which fits a logistic regression of
  #   these parameters on the proportion of occupied neighbouring cells
  #   (in a queen's or 2nd order neighbourhood) during the previous time step
  # (5) Additional detection heterogeneity can be introduced
  #   at the site- or the individual survey level, with the possibility of a
  #   temporal trend in this heterogeneity. For instance, an annual trend in
  #   detection heterogeneity at the site or the survey level is specified by
  #   the value in the first and the last year.
  #   Hence, trend.sd.site = c(0, 1) will result in a linear trend in
  #   the magnitude of site heterogeneity in detection from 0 in the
  #   first year to 1 in the last year.
  #
  #
  # Function arguments:
  # -------------------
  #
  # *** Design of study and basic 'magnitude' of parameters ***
  # side – side length of square simulation area. Therefore,
  #     the number of sites, or cells, M = side^2
  # nsurveys – Number of replicate surveys within a 'season', year or primary period
  # nyears – Number of years (or 'seasons')
  # mean.psi1 – intercept of occupancy probability in year 1
  # range.phi and range.gamma – bounds of uniform distribution from which
  #   annual intercepts for persistence (phi) and colonisation (gamma)
  #   are drawn
  # range.p – same for detection probability p
  #
  #
  # *** Covariates ***
  # beta.Xpsi1: coefficient of a site covariate in psi1
  # beta.Xphi: coefficient of a site/year covariate in phi
  # beta.Xgamma: coefficient of a site/year covariate in gamma
  # beta.Xp: coefficient of a site/year/rep covariate in p
  #
  #
  # *** Parameters governing the spatial correlations ***
  # theta.XAC: 'range parameter' of a covariate with exponential
  #   spatial correlation (i.e., a Gaussian random field is used as an
  #   environmental covariate). NOTE: if you want to set to zero the effects
  #   of this spatially autocorrelated variable, you CANNOT
  #   set theta.XAC=0 because this breaks the function,
  #   nor can you simply choose a very small value.
  #   Instead you MUST set the elements of coefficients vector beta.XAC
  #   to zero.
  # beta.XAC: vector of coefficients of that field for the 4 model params:
  #   psi1, phi, gamma, and p (in that order)
  # beta.Xautolog – vector of coefficients of autologistic covariate
  #   in the following order: persistence (phi), colonization (gamma).
  #   Autocovariate is computed at every season as the proportion of
  #   occupied cells in a queen's neighbourhood around each cell.
  #
  #
  # *** Detection heterogeneity ***
  # trend.sd.site: range of year-specific values of SD of Gaussian
  #   random site effects in p: c(1,1) specifies constant value of 1
  #   for all years, while c(0,1) specifies linear increase over the years
  #   from 0 to 1.
  # trend.sd.survey: range of year-specific values of standard deviation
  #   of Gaussian random survey effects in p: specification as
  #   for trend.sd.site
  #
  # *** Graphics control and other ***
  # seed – allows to 'fix' the simulation such that it becomes reproducible
  # ask.plot – if TRUE permits to browse through plots (otherwise if FALSE)

  if(FALSE) {x <- NULL; rm(x)}  # Stops R CMD check choking on 'curve'.

  # Checks and fixes for input data -----------------------------
  side <- round(side[1])
  nyears <- round(nyears[1])
  stopifnotGreaterthan(nyears, 1)
  nsurveys <- round(nsurveys[1])
  stopifnotProbability(mean.psi1)
  stopifnotProbability(range.phi) # bounds
  stopifnotProbability(range.gamma) # bounds
  stopifnotProbability(range.p) # bounds
  stopifNegative(theta.XAC, allowZero=FALSE)
  stopifnotLength(beta.XAC, 4)
  stopifnotLength(beta.Xautolog, 2)
  stopifnotLength(trend.sd.site, 2) # trend
  stopifNegative(trend.sd.site)
  stopifnotLength(trend.sd.survey, 2) # trend
  stopifNegative(trend.sd.survey)
  # ----------------------------------------------------------------
  # Restore graphical settings on exit -----------------------------
  if(show.plots) {
    oldpar <- par("mfrow", "mar", "cex.main", "cex.lab", "cex.axis")
    oldAsk <- devAskNewPage(ask = ask.plot && dev.interactive(orNone=TRUE))
    on.exit({par(oldpar); devAskNewPage(oldAsk)})
  }
  # ----------------------------------------------------------------

  # Create grid
  xcoord <- 1:side
  ycoord <- 1:side
  grid <- as.matrix(expand.grid(x=xcoord, y=ycoord))
  M <- side^2             # Total number of cells or sites

  # Compute adjacency matrix for grid
  neigh <- spdep::dnearneigh(as.matrix(grid), d1 = 0, d2 = sqrt(2) * 1 + 0.01)
  winnb <- spdep::nb2WB(neigh)        # Function to get CAR ingredients for BUGS
  nneigh <- winnb$num          # number of neighbours
  amatrix <- spdep::nb2mat(neigh)
  amatrix[amatrix > 0] <- 1    # Neighbours get a 1, non-neighbours a 0

  # Set up arrays needed
  site <- 1:M                                      # Sites
  year <- 1:nyears                                 # Years
  prob <- array(dim = c(side, side))               # p matrix
  psi <- muZ <- z <- array(dim = c(side, side, nyears))      # Occupancy, occurrence
  phi <- gamma <- array(NA, dim = c(side, side, (nyears-1))) # Survival, colonisation
  Xauto <- array(NA, dim = c(side, side, nyears))            # Autocovariate
  y <- p <- array(NA, dim = c(side, side, nsurveys, nyears)) # Det. histories and p

  # Create values of 1 spatially autocorrelated covariate XAC
  # Generate correlated random variables in a square
  RandomFields::RFoptions(seed=seed.XAC)     # Default NA; 88 gives cool pattern
  XAC <- matrix(RandomFields::RFsimulate(RandomFields::RMexp(var = 1, scale = theta.XAC),
      x=xcoord, y=ycoord, grid=TRUE)@data$variable1,
      ncol = side, byrow = TRUE)  # variance 1
  if(!is.na(seed.XAC))
    RandomFields::RFoptions(seed=NA)

  set.seed(seed=seed)  # Default NULL; do this AFTER RFsimulate

  # Create four spatially unstructured covariates
  # Site covariate for psi1
  Xpsi1 <- matrix(runif(M, -2, 2), ncol = side)

  # Yearly-site covariates for phi and gamma
  Xphi <- Xgamma <- array(runif(M*nyears, -2, 2), dim = c(side, side, nyears))

  # Observational covariate for p
  Xp <- array(runif(M*nsurveys*nyears,-2,2), dim=c(side, side,nsurveys,nyears))

  # Draw values of baseline levels of the main parameters
  # (i.e., draw year effects if any)
  mean.phi <- runif(n = nyears-1, min = min(range.phi), max = max(range.phi))
  mean.gamma <- runif(n = nyears-1, min = min(range.gamma), max = max(range.gamma))
  mean.p <- runif(n = nyears, min = min(range.p), max = max(range.p))



  # (a) Simulate state process parameters: initial state (first year)
  psi[,,1] <- plogis(qlogis(mean.psi1) + beta.Xpsi1 * Xpsi1 +
     beta.XAC[1] * XAC) # psi1

  # (b) Simulate state in first year
  z[,,1] <- rbinom(M, 1, psi[,,1])  # Initial occurrence state

  # Compute value of autocovariate after first year = proportion of neighbours occupied
  # first vectorize and then put into matrix again
  Xautovec <- amatrix %*% c(z[,,1])
  Xauto[,,1] <- matrix(Xautovec/nneigh, ncol = side) # Put back in matrix by column again
  
  # Do the pre-loop plots
  # ---------------------
  if(show.plots) {
    tryPlot <- try( {
      # Plot effects of autocovariate on (year-specific) phi and gamma
      par(mfrow = c(1, 2))
      curve(plogis(qlogis(mean.phi[1]) + beta.Xautolog[1] * x), 0, 1, 
          main = "Persistence: \nphi ~ Year + Autocovariate", xlab = "Autocov. (prop. occupied neighb.)",
          ylab = "phi", ylim = c(0,1), frame = FALSE)
      for(k in 2:(nyears-1)){
         curve(plogis(qlogis(mean.phi[k])+beta.Xautolog[1]*x),0,1,add=TRUE)
      }

      curve(plogis(qlogis(mean.gamma[1]) + beta.Xautolog[2] * x), 0, 1, 
          main = "Colonization: \ngamma ~ Year + Autocovariate",
          xlab = "Autocovariate (prop. occupied neighb.)", ylab = "gamma", ylim = c(0,1), frame = FALSE)
      for(k in 2:(nyears-1)){
         curve(plogis(qlogis(mean.gamma[k])+beta.Xautolog[2]*x),0,1,add=TRUE)
      }

      # Simulate true system dynamics
      par(mfrow = c(2,2), mar = c(5,4,5,2), cex.main = 1.3, cex.lab = 1.5, cex.axis = 1.2)

      # Plot random field covariate XAC
      # rows are in x, columns in y direction
      image(1:side, 1:side, XAC, col=topo.colors(100), 
          main = paste("Gaussian random field XAC with \n neg. exponential correlation (range =", theta.XAC, ")"), 
          xlab = 'x', ylab = 'y')

      image(1:side, 1:side, psi[,,1], col=topo.colors(100),
          main = paste("Initial occupancy probability"), xlab = 'x', ylab = 'y')

      image(1:side, 1:side, z[,,1], col=c("white", "black"), 
          main = paste("Initial presence/absence (true system state z):\n black = occupied, white = unoccupied"),
          xlab = 'x', ylab = 'y')
      abline(h = 0:side+0.5, v = 0:side+0.5, col = "lightgrey")

      image(1:side, 1:side, Xauto[,,1], col=topo.colors(100), 
          main = "Autocovariate between year 1 and year 2", xlab = 'x', ylab = 'y')
    }, silent = TRUE)
    if(inherits(tryPlot, "try-error")) {
      show.plots <- FALSE
      tryPlotError(tryPlot)
    }
  }

  # (c) Simulate state process parameters: time steps 2:nyears
  for(k in 2:nyears){
    par(mfrow = c(2,2), mar = c(5,4,5,2), cex.main = 1.3, cex.lab = 1.5, cex.axis = 1.2)
    if(verbose)
      cat(paste("** Year", k, "**\n"))
    # Compute colonisation and extinction parameters and plot
    phi[,,k-1] <- plogis(qlogis(mean.phi[k-1]) + beta.Xphi * Xphi[,,k-1] + beta.XAC[2] * XAC +
      beta.Xautolog[1] * Xauto[,,k-1])
    gamma[,,k-1] <- plogis(qlogis(mean.gamma[k-1]) + beta.Xgamma * Xgamma[,,k-1] + beta.XAC[3] * XAC +
      beta.Xautolog[2] * Xauto[,,k-1])

    # Compute latent states and plot
    muZ[,,k] <- z[,,k-1]*phi[,,k-1] + (1-z[,,k-1])*gamma[,,k-1]
    z[,,k] <- rbinom(M, 1, muZ[,,k])

    # Compute autocovariate and plot
    Xautovec <- amatrix %*% c(z[,,k])
    Xauto[,,k] <- matrix(Xautovec/nneigh, ncol = side) # re-assemble by column
    
    # Do the in-loop plots
    # --------------------
    if(show.plots) {
      tryPlot <- try( {
        image(1:side, 1:side, phi[,,k-1], col=topo.colors(100),
              main = paste("Persistence between year", k-1, "and year", k), xlab = 'x', ylab = 'y')
        image(1:side, 1:side, gamma[,,k-1], col=topo.colors(100),
              main = paste("Colonization between year", k-1, "and year", k), xlab = 'x', ylab = 'y')
        image(1:side, 1:side, z[,,k], col=c("white", "black"),
            main = paste('Presence/absence (z) in year', k, ':\n black = occupied, white = unoccupied'),
            xlab = 'x', ylab = 'y')
        abline(h = 0:side+0.5, v = 0:side+0.5, col = "lightgrey")
        image(1:side, 1:side, Xauto[,,k], col=topo.colors(100),
            main = paste("Autocovariate between year", k, "and year", k+1), xlab = 'x', ylab = 'y')
      }, silent = TRUE)
      if(inherits(tryPlot, "try-error")) {
        show.plots <- FALSE
        tryPlotError(tryPlot)
      }
    }

  }

  # (d) Observation process parameters
  # First choose values of annual SD of p random effects
  sd.site <- seq(from = trend.sd.site[1], to = trend.sd.site[2],
     length.out = nyears)
  sd.survey <- seq(from = trend.sd.survey[1], to = trend.sd.survey[2],
     length.out = nyears)
  for(k in 1:nyears){
    # Site random effects
    eps1 <- matrix(rnorm(n = M, mean = 0, sd = sd.site[k]), ncol = side)
    # Survey random eff.
    eps2 <- rnorm(n = nsurveys, mean = 0, sd = sd.survey[k])
    for(j in 1:nsurveys){
      p[,,j,k] <- plogis(qlogis(mean.p[k]) + beta.Xp * Xp[,,j,k] + beta.XAC[4] * XAC + eps1[,] + eps2[j])
    }
  }

  # Simulate actual observation process (also updating entire grid in one go)
  for(k in 1:nyears){                 # Loop over years
    for(j in 1:nsurveys){               # Loop over replicates
      prob <- z[,,k] * p[,,j,k]  # zero out p for unoccupied sites
      y[,,j,k] <- rbinom(M, 1, prob)
    # image(1:side, 1:side, y[,,j,k], main = paste("Year", k, "and rep", j))
        # Look at clumped pattern in y
    }
  }

  # Derived quantities
  # Compute annual population occupancy
  for (k in 2:nyears){
     psi[,,k] <- psi[,,k-1]*phi[,,k-1] + (1-psi[,,k-1])*gamma[,,k-1]
  }
  mean.psi <- apply(psi, 3, mean)              # Average psi over sites

  # Compute true and observed number of occupied sites
  zobs <- apply(y, c(1,2,4), max)
  nocc <- apply(z, 3, sum)
  nocc.obs <- apply(zobs, 3, sum)

  # Do the post-loop plots
  # ----------------------
  psi.app <- apply(apply(y, c(1,2,4), max), 3, mean)
  if(show.plots) {
    tryPlot <- try( {
      # (4) More plots comparing true states and observations
      # Plot realised and apparent occupancy
      par(mfrow = c(1,1))
      plot(year, apply(z, 3, mean), type = "l", xlab = "Year", ylab = "Occupancy or Detection prob.", col = "red", xlim = c(0,nyears+1), ylim = c(0,1), lwd = 2, lty = 1, frame.plot = FALSE, las = 1)
      lines(year, mean.p, type = "l", col = "red", lwd = 2, lty = 2)
      lines(year, psi.app, type = "l", col = "black", lwd = 2)
      text(0.8*nyears, 0.1, labels = "red solid - true occupancy prob.\n red dashed - detection prob.\n black - observed proportion occupied", cex = 1)

      # Plots comparing true and observed latent states
      par(mfrow = c(2,2), mar = c(5,4,5,2), cex.main = 1.3, cex.lab = 1.5, cex.axis = 1.2)
      for(k in 1:nyears){
        image(1:side, 1:side, z[,,k], col=c("white", "black"), main = paste('Presence/absence (z) in year', k), xlab = 'x', ylab = 'y')
        abline(h = 0:side+0.5, v = 0:side+0.5, col = "lightgrey")
        image(1:side, 1:side, zobs[,,k], col=c("white", "black"), main = paste('Ever_detected (zobs) in year', k), xlab = 'x', ylab = 'y')
        abline(h = 0:side+0.5, v = 0:side+0.5, col = "lightgrey")
      }
    }, silent = TRUE)
    if(inherits(tryPlot, "try-error"))
      tryPlotError(tryPlot)
  }
  
  # Compute values of naive autocovariate (observed prop.
  #     of occupied neighbouring cells)
  Xautoobs <- array(NA, dim = dim(zobs))
  for(k in 1:nyears){            # Loop over years
    for(i1 in 1:side){      # Loop over one side (hopefully X)
      for(i2 in 1:side){    # Loop over other side (hopefully Y)
        i1.start <- max(1,(i1-1))
        i1.end <- min(side,(i1+1))
        i2.start <- max(1,(i2-1))
        i2.end <- min(side,(i2+1))
        Xautoobs[i1,i2,k] <- (sum(zobs[i1.start:i1.end,i2.start:i2.end,k])- zobs[i1,i2,k])  / (length(zobs[i1.start:i1.end,i2.start:i2.end,k]) - 1)
      }
    }
  }

  out <- list(
    # ----------------- values input -----------------------
    side=side, nyears=nyears, nsurveys=nsurveys, mean.psi1=mean.psi1, beta.Xpsi1=beta.Xpsi1,
    range.phi=range.phi, beta.Xphi=beta.Xphi,
    range.gamma=range.gamma, beta.Xgamma=beta.Xgamma,
    range.p=range.p, beta.Xp=beta.Xp,
    theta.XAC=theta.XAC, beta.XAC= beta.XAC, beta.Xautolog=beta.Xautolog,
    trend.sd.site=trend.sd.site, trend.sd.survey=trend.sd.survey,
    seed=seed, seed.XAC = seed.XAC,
    # ----------------- values generated --------------------
    M=M,                # total number of pixels in the study area
    grid=grid,          # 2-column matrix, the x and y coordinates of the pixels
    amatrix = amatrix,  # MxM matrix, [i,j] = 1 if cell i and j are neighbours, 0 otherwise
    Xpsi1 = Xpsi1,      # side x side matrix, value of covariate affecting initial occupancy (psi1)
    Xphi = Xphi,        # side x side x nyears array, value of covariate affecting persistence (phi)
    Xgamma = Xgamma,    # side x side x nyears array, value of covariate affecting colonisation (gamma)
    Xp = Xp,            # side x side x nsurveys x nyears array, value of covariate affecting detection (p)
    XAC=XAC,            # side x side matrix, the spatially correlated covariate
    Xauto = Xauto,      # side x side x nyears array, autocovariate, proportion of neighbouring cells occupied
    Xautoobs = Xautoobs,# side x side x nyears array, observed autocovariate, proportion of neighbouring cells where species detected
    sd.site=sd.site,    # vector nyears, year-specific values of SD of Gaussian random site effects in p
    sd.survey=sd.survey,# vector nyears, year-specific values of SD of Gaussian random survey effects in p
    mean.phi=mean.phi,  # vector nyears-1, year-specific intercept of persistence on probability scale
    mean.gamma=mean.gamma,# vector nyears-1, year-specific intercept of colonisation on probability scale
    mean.p=mean.p,      # vector nyears, year-specific intercept of detection probability on probability scale
    psi=psi,            # side x side x nyears array, probability of occupancy of cell
    mean.psi=mean.psi,  # vector nyears, mean occupancy over all cells
    psi.app=psi.app,    # vector nyears, apparent occupancy, proportion of cells where species detected
    z=z,                # side x side x nyears array, true occupancy status of each cell in each year (1 if occupied)
    zobs=zobs,          # side x side x nyears array, observed occupancy status of each cell in each year (1 if detected)
    nocc = nocc,        # vector nyears, the true number of cells occupied each year
    nocc.obs = nocc.obs,# vector nyears, the number of cells where species detected each year
    phi=phi,            # side x side x nyears-1 array, probability of persistence in each interval between years
    gamma=gamma,        # side x side x nyears-1 array, probability of colonisation in each interval between years
    p=p,                # side x side x nsurveys x nyears array, probability of detection
    y = y)             # side x side x nsurveys x nyears array, detection history, 1 if species detected.

  # Add an unmarked data frame object
  out$umf <- conv2UM(out)

  return(out)
} # ------------------------------------------------------------------

# Helper function that turns the simulated data into an unmarked data frame
conv2UM <- function(d){
# Function creates an unmarked data frame for model-fitting function colext
#   (= dynamic occupancy model) using output of the data simulation function.
#
# The function works for up to 9999 years = 'seasons' or 'primary periods'.
#
# Marc Kery, 12 Dec 2014; mangled by Mike Meredith, 17 May 2019
#
# Function arguments:
# d: output object of the data sim function

# Row names to be used for everything:
names <- paste(d$grid[,1], d$grid[,2], sep = '.')

yy <- matrix(d$y, d$M, d$nsurveys * d$nyears)
Xp <- matrix(d$Xp, d$M, d$nsurveys * d$nyears)
rownames(yy) <- rownames(Xp) <- names

siteCovs <- data.frame(Xpsi1 = c(d$Xpsi1), XAC = c(d$XAC))
rownames(siteCovs) <- names

if(d$nyears < 100) {
  yrChar <- sprintf("%02i", 1:d$nyears)
} else {
  yrChar <- sprintf("%04i", 1:d$nyears)
}

yearlySiteCovs <- list(
    year = matrix(yrChar, d$M, d$nyears, byrow=TRUE, dimnames=list(names, NULL)),
    Xphi = matrix(d$Xphi, d$M, d$nyears, dimnames=list(names, NULL)),
    Xgamma = matrix(d$Xgamma, d$M, d$nyears, dimnames=list(names, NULL)),
    Xauto = matrix(d$Xauto, d$M, d$nyears, dimnames=list(names, NULL)),
    Xautoobs = matrix(d$Xautoobs, d$M, d$nyears, dimnames=list(names, NULL)) )

# Create and return unmarked data frame for the colext function
return(unmarked::unmarkedMultFrame(y=yy, siteCovs=siteCovs, yearlySiteCovs=yearlySiteCovs,
  obsCovs = list(Xp = Xp), numPrimary = d$nyears))
} # --------------------------------------
