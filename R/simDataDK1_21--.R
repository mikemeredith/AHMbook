
# A function to simulate data under a (thinned) inhomogenous Poisson point process (IPP)
# Code adapted from Dorazio (GEB, 2014) and Koshkina et al. (MEE, 2017)
#  by Marc then decorticated by Mike.

# Like the original D-K code, this version only allows one individual per pixel.

# AHM2 chapter 10

# Helper function to generate bivariate normal covariate surfaces; not exported:
# getCovSurface <- function(mWt=c(0.75, 0.4),sWt=c(0.25, 0.5), rho=0.5,
    # xmin=-1, xmax=1, ymin=-1, ymax=1, loc) {
  # mu.x <- xmin + mWt[1]*(xmax-xmin)
  # mu.y <- ymin + mWt[2]*(ymax-ymin)
  # sigma.x <- sWt[1]*abs(xmax-xmin)
  # sigma.y <- sWt[2]*abs(ymax-ymin)
  # rho1.xy <- 0.5
  # mu <- c(mu.x, mu.y)
  # Sigma <- matrix(c(sigma.x^2,
      # rep(rho*sigma.x*sigma.y, 2),
      # sigma.y^2), ncol=2)
  # mvtnorm::dmvnorm(loc, mean=mu, sigma=Sigma)
# }  # ----------------------------------------



# ---------------- Start of function definition --------------------
simDataDK1 <- function(sqrt.npix = 100, alpha = c(-1,-1), beta = c(6,0.5),
  drop.out.prop.pb = 0.7, quadrat.size = 4, gamma = c(0,-1.5),
  nquadrats = 250, nsurveys = 3, show.plot = TRUE){
  #
  # Function generates data for use with the integrated model described by Dorazio (GEB, 2014)
  # The function is based on the code written by Dorazio and adapted by Koshkina et al. (MEE, 2017)
  #
  # A Poisson point pattern (PPP) with intensity a function of a covariate X and intercept and
  # coefficient beta is simulated on a discrete (pixel-based) approximation of a continuous landscape
  # This PPP is then thinned with a pixel-wise thinning probability and with a
  # landscape-wise drop.out.prob.pb to produce a first data set of presence-only kind.
  # A second data set is simulated by imagining replicated pount counts in a total of nquadrats quadrats among
  # all square partitions of the landscape by quadrat.size.


  # Simulate data under a (thinned) Poisson point pattern (PPP).
  # sqrt.npix: number of pixels along each side of square state space (the 'landscape')
  #    number of pixels is then sqrt.npix^2
  # intensity of IPP: log(lambda) = beta0 + beta1 * covariate X
  # sampling detection bias in presence-only observations of IPP is modelled as:
  #               logit(b) = alpha0 + alpha1 * covariate W
  # quadrat.size: length of the side of quadrats in *pixel* units by which PPP is summarized for conducting replicate counts or site-occ surveys
  # detection probability for the counts is governed by gamma0 and gamma1 (on the cell-averaged values of W)
  # Drop-out proportion is proportion of PO points at the end that are discarded (perhaps because these sites are not visited at all ?)
  # nsurveys: number of replicated surveys in the count survey


  # This is a 100x smaller version of the simulation field compared to Koshkina et al.
  # The grid is only 100 x 100 = 10k (instead of 1 M)

  # For the purpose of getting count surveys, we aggregate the original landscape by quadrats of size 16, which will yield 10k / 16 = 625 cells of 4 by 4 units side length

  # Of these, we will keep a random sample of 'nquadrats' quadrats for the count surveys.

  # Note that the parameters beta must be chosen such to avoid a too high 'filling' of the discrete approximation of the entire field B
  # ------------------------------------------------------

  # -------------- Check and fix input -----------------------
  sqrt.npix <- round(sqrt.npix[1])
  stopifnotLength(alpha, 2)
  stopifnotLength(beta, 2)
  stopifnotProbability(drop.out.prop.pb)
  quadrat.size <- round(quadrat.size[1])
  if(sqrt.npix %% quadrat.size != 0)
    stop("sqrt.npix / quadrat.size must return an integer.", call.=FALSE)
  stopifnotLength(gamma, 2)
  nquadrats <- round(nquadrats[1])
  if(nquadrats > (sqrt.npix / quadrat.size)^2)
    stop("Number of quadrats to sample exceeds number of quadrats in the landscape.", call. = FALSE)
  nsurveys <- round(nsurveys[1])
  # ------------------------------------------------------------


  # Define landscape as a rectangular region S, with x and y ranging
  #  from -1 to +1, hence:
  s.area <- 4

  # Approximate the landscape by many small pixels in a raster object
  s <- temp <- raster(ncol=sqrt.npix, nrow=sqrt.npix, xmn=-1, xmx=1, ymn=-1, ymx=1, crs=NULL)
  # s will become a Raster Stack, temp is a template to create new layers
  s.loc <- xyFromCell(temp, 1:ncell(s))    # Coordinates of every pixel in S

  # PART A. Ecological process - where are the animals?
  # ---------------------------------------------------
  # Compute covariate X as the sum of two MVN random variables
  xcov <- 0.4 * getCovSurface(mWt=c(0.75,0.4), sWt=c(0.25, 0.5),rho=0.5, loc=s.loc) +
          0.6 * getCovSurface(mWt=c(0.15,0.8), sWt=c(0.5, 0.25),rho=-0.4, loc=s.loc)
  xcov <- standardize(xcov)    # Standardize covariate x

  # Fill covariate x into the raster s
  values(s) <- xcov
  names(s) <- 'x'

  # Calculate lambda as a function of X and add to Raster Stack
  X <- cbind(1, xcov)
  values(temp) <- exp(X %*% beta)   # calculate lambda
  names(temp) <- 'lambda'
  s <- addLayer(s, temp)  # add to raster stack

  # We use rejection sampling to get random draws from the IPP. The proposal
  #  distribution is uniform, ie, HPP:
  # How 'high' should the proposal density be? It must exceed the IPP everywhere.
  ( maxlambda <- max(values(s)[,'lambda']) ) # ~ 900 per unit area
  (N.hpp <- suppressWarnings(rpois(1, maxlambda*s.area)))     # Number we need to draw
  if(is.na(N.hpp) || N.hpp >= ncell(s))
    stop("The 'beta' settings result in intensities that are too high\nin the most intense region (more animals than pixels.)", call.=FALSE)
  # Draw from the proposal distribution
  ind.hpp <- sample(1:ncell(s), size = N.hpp, replace = FALSE)   #  sampling w/o replacement ensures only 1 individual per pixel
  # reject draws depending on lambda (and hence X) to get the IPP draws
  lambda.hpp <- values(s)[,'lambda'][ind.hpp]  # intensities at those pixels
  ind.ipp <- rbinom(N.hpp, 1, lambda.hpp/maxlambda)   # use intensity lambda to determine whether to accept or reject.
  (N.ipp <- sum(ind.ipp))               # about 1800 individuals in IPP
  pixel.id.ipp <- ind.hpp[ind.ipp == 1] # Gives id of every pixel in landscape that has an individual, a vector of numbers, length N.ipp = total population

  # PART B. The detection-only observation model
  # --------------------------------------------
  # Create covariate W which affects which animals are detected (ie, where people look)
  wcov <- getCovSurface(mWt=c(0.25,0.65), sWt=c(0.25, 0.5),rho=0.1, loc=s.loc)
  wcov <- standardize(wcov)            # Standardize covariate w
  # Add to the raster stack
  values(temp) <- wcov
  names(temp) <- 'w'
  s <- addLayer(s, temp)

  # Compute value of thinning parameter b (= probability of detection)
  #  for each cell as a function of alpha and covariate W
  W <- cbind(1, wcov)   # Design matrix W for thinning
  values(temp) <- plogis(W %*% alpha) # thinning coefs
  names(temp) <- 'pTrue'
  s <- addLayer(s, temp)

  # ... simulate presence-only data (= detections of individuals as a thinned point process)
  pTrue.ipp <- values(s)[,'pTrue'][pixel.id.ipp]
  y.ipp1 <- rbinom(N.ipp, size=1, prob=pTrue.ipp) # 1 = detected, 0 = not detected
  # A 1/0 vector, length N.ipp
  pixel.id.det1 <- pixel.id.ipp[y.ipp1 == 1]
  # length(pixel.id.det1) # ~ 550. This is too many, better drop a bunch

  drop.out <- runif(length(pixel.id.det1), 0, 1) < drop.out.prop.pb # T/F vector
  pixel.id.det <- pixel.id.det1[!drop.out]
  # length(pixel.id.det) # ~ 180, ok
  y.point <- numeric(sqrt.npix^2)
  y.point[pixel.id.det] <- 1


  #### Part C: simulate replicate count data
  ## ---------------------------------------
  # Get the info on animal location into our Raster Stack as a layer:
  spop <- dropLayer(s, c('lambda','pTrue')) # clean up, just keep covars
  z <- rep(0, ncell(spop))
  z[pixel.id.ipp] <- 1
  values(temp) <- z
  names(temp) <- 'presence'  # 1 if animal in pixel, 0 otherwise (max 1 animal per pixel)
  spop <- addLayer(spop, temp)

  # Form quadrats of side quadrat.size (default 4 -> area 16 -> 625 quadrats)
  quadfact <- c(quadrat.size, quadrat.size)
  squad <- raster::aggregate(spop, fact=quadfact, fun=mean)
  # mean is ok for x and w, but for 'presence' we need the sum
  abund <- raster::aggregate(raster::subset(spop, 'presence'), fact=quadfact, fun=sum)
  names(abund) <- 'N'
  squad <- addLayer(squad, abund)
  squad <- dropLayer(squad, 'presence') # clean up, N/16 not useful

  # Simulate replicate counts at every quadrat (aka "site")
  nsite <- ncell(squad)          # number of sites/quadrats in count design
  N <- values(squad)[,'N']       # Extract latent abundance at each site

  # Compute values of detection probability for each 16-cell pixel in count survey
  # (Here we use the same covar, w, as the detection-only model, as do Koshkina et al,
  #   but a different covar could be generated.)
  pcount <- plogis(gamma[1] + gamma[2] * values(squad)[,'w'])

  # Do the nsurveys surveys
  counts <- array(NA, dim = c(nsite, nsurveys))
  for(j in 1:nsurveys){
    counts[,j] <- rbinom(nsite, N, pcount)
  }

  fullCountData <- cbind(quadID=1:nsite, values(squad), counts)

  # Draw random sample of 'nquadrats' quadrats from the total number, nsite
  selQuad <- sort(sample(1:nsite, nquadrats, replace = FALSE))
  countData <- fullCountData[selQuad,]

  # Output (visual): shown only when show.plot = TRUE
  if(show.plot) {
    oldpar <- par(mfrow=1:2, mar = c(2,1,6,3))
    oldAsk <- devAskNewPage(ask = dev.interactive(orNone = TRUE))
    on.exit({par(oldpar) ; devAskNewPage(oldAsk)})

    # Fig.1
    loc.ipp <- s.loc[pixel.id.ipp, ]
    raster::plot(raster::subset(s, 'x'), axes = FALSE, box = FALSE, asp=1,
      main = paste("Inhomogenous Poisson point process:\nIntensity  covariate 'x' and\nlocations of", N.ipp, "individuals"))
    points(loc.ipp, pch = 16, cex = 0.5)    # location of the individuals

    loc.det <- s.loc[pixel.id.det, ]
    N.det <- length(pixel.id.det)
    raster::plot(raster::subset(s, 'w'), axes = FALSE, box = FALSE,  asp=1,
      main = paste("Presence-only observations:\nDetection bias covariate 'w' and\nlocations of", N.det, "individuals detected"))
    points(loc.det, pch = 16, cex = 0.5)    # location of the individuals

    # Fig.2
    par(mfrow=c(2,2), mar = c(1,1,5,3))
    raster::plot(raster::subset(squad, 'x'), axes = FALSE, box = FALSE,  asp=1,
      main = "Mean intensity covariate 'x'\nfor each quadrat")
    raster::plot(raster::subset(squad, 'w'), axes = FALSE, box = FALSE,  asp=1,
    main = "Mean detection covariate 'w'\nfor each quadrat")
    raster::plot(raster::subset(squad, 'N'), axes = FALSE, box = FALSE,  asp=1,
      main = "True abundance 'N'\nfor each quadrat" )
    mnc <- rowMeans(counts)
    mnc[-selQuad] <- NA
    cnt <- raster::subset(squad, 'N')
    values(cnt) <- mnc
    raster::plot(cnt, colNA='darkgrey',axes = FALSE, box = FALSE, asp=1,
      main = "Mean counts for \neach quadrat surveyed,\ngrey if unsurveyed")
  }    # end show.plot

  # Output (numeric)
  return(list(
    # ---------------- input arguments ------------------
    sqrt.npix = sqrt.npix, alpha = alpha, beta = beta, gamma = gamma,
    drop.out.prop.pb = drop.out.prop.pb, quadrat.size = quadrat.size,
    nquadrats = nquadrats, nsurveys =nsurveys,
    # ----------------- values generated -------------------
    npix = sqrt.npix^2,            # Number of pixels in the landscape
    s.area = s.area,               # Area of the landscape, 4
    s.loc = s.loc,                 # Coordinates of every pixel in the landscape
    xcov = xcov,                   # 'x' (intensity) covariate
    wcov = wcov,                   # 'w' (detection) covariate
    N.ipp = N.ipp,                 # True number of individuals in the landscape
    pixel.id.ipp = pixel.id.ipp,   # Pixel ID for each individual in the population
    loc.ipp = s.loc[pixel.id.ipp, ],
        # Coordinates for each individual in the population
    pTrue.ipp = pTrue.ipp,         # Probability of detection for each individual
    pixel.id.det = pixel.id.det,   # Pixel ID for each individual detected
    N.det = length(pixel.id.det),  # Number of detections
    loc.det = s.loc[pixel.id.det, ], # Coordinates for each individual detected
    pcount = pcount,               # Probability of detection in each quadrat
    fullCountData = fullCountData,
      # matrix with rows for each quadrat, columns for ID, x and w coords,
      #   true N, and 3 replicate counts
    countData = countData,         # as above, but rows for quadrats sampled only
    s = s,                         # Raster Layer object with ... for all pixels
    squad = squad))                # Raster Layer object with ... for quadrats
} # ---------------- End of function definition --------------------

