
# Function to generate the Habitat object used in section 11.10
# Density Surface Modeling

# X : a 2-column matrix with coordinates of _regularly_spaced_ points along the transect line
#   (X <- regpoints@coords in book)
# Ntotal : true number of individuals in the study area
# sigma.move = 0 :  not used! I took it out ####
# sigma : scale parameter for the half-normal detection function
# beta1 : coefficient for the relationship between density and the habitat covariate
# nSurveys : the number of surveys to simulate
# xlim, ylim : the extent of the (rectangular) study area

simDSM <- function(X, Ntotal = 400, sigma = 0.65, beta1 = 1.0,
    nsurveys = 2, xlim = c(-0.5, 3.5), ylim = c(-0.5, 4.5), show.plots = TRUE) {

  # Create coordinates rasterized transect
  delta <- 0.2               # 2D bin width
  gry <- seq(ylim[1] + delta/2, ylim[2] - delta/2, delta)
  ny <- length(gry)
  grx <- seq(xlim[1] + delta/2, xlim[2] - delta/2, delta)
  nx <- length(grx)
  grx <- rep(grx, ny)
  gry <- rev(sort(rep(gry,nx)))
  gr <- cbind(grx, gry)
  nPix <- nrow(gr)
  # Create spatially correlated covariate x and plot it
  V <- exp(-e2dist(gr, gr)/1)
  x <- t(chol(V)) %*% rnorm(nrow(gr))  ### x = habitat covariate, stochastic

  # Simulate activity centre locations
  probs <- exp(beta1*x)/(sum(exp(beta1*x)))
  # Activity centers selected based on habitat
  s.pix.id <- sample(1:nPix, Ntotal, prob = probs, replace=TRUE)  ### stochastic - sample ###
  N <- tabulate(s.pix.id, nbins = nPix)
  # Uniformly distributed within their pixel
  sx <- runif(Ntotal, gr[s.pix.id,1] - delta/2, gr[s.pix.id,1] + delta/2)  ### stochastic
  sy <- runif(Ntotal, gr[s.pix.id,2] - delta/2, gr[s.pix.id,2] + delta/2)
  U <- cbind(sx,sy)
  ### new simulation - observation process

  # Compute p for each activity center for each point along the line
  parr <- numeric(Ntotal) # same for all surveys
  y2d <- array(0, c(Ntotal, nsurveys))  # 0/1 detection matrix
  for (i in 1:Ntotal) {
    dvec <- min( sqrt((sx[i] - X[, 1])^2 +(sy[i] - X[, 2])^2) )
    loghaz <- -(1/(2*sigma*sigma)) * dvec * dvec
    parr[i] <- exp(loghaz)
    y2d[i, ] <- rbinom(nsurveys, 1, parr[i] )  ### stochastic
  }

  # Find which individuals are captured at least once, REMOVE the rest
  cap <- apply(y2d, 1, sum) > 0   # TRUE  if captured at least once
  nind <- sum(cap)                # number captured at least once
  y2d <- y2d[cap, , drop=FALSE]
  Ucap <- U[cap, ]  # matrix with AC coords
  gid <- s.pix.id[cap]      # pixel IDs for ACs

  # generate pixel matrix with nsurvey columns, with NA if not captured
  pixel <- matrix(gid, nrow=nind, ncol=nsurveys)
  pixel[y2d == 0] <- NA
  # for back compatibility we need to reverse the order and add colnames
  pixel <- pixel[nind:1, , drop=FALSE]
  colnames(pixel) <- paste0("gid", 1:nsurveys)
  # End of data preparation

  # Plot with activity centers and linking lines # Fig. 11-15
  if(show.plots) {
    oldpar <- par(mar = c(3,3,3,6)); on.exit(par(oldpar))
    tryPlot <- try( {
      image(r <- rasterFromXYZ(cbind(gr,x)), col = topo.colors(10))
      image_scale(x, col = topo.colors(10))
      lines(X, col = "black", pch = 20, lwd = 3)
      points(sx, sy, pch = 16, col = "black", lwd=1 )
      # plot observed locations (i.e. detected individuals)
      points(Ucap, pch = 20, col = "red")
      # Add lines from detected ACs to nearest point on transect
      dd <- e2dist(X, Ucap)
      closest <- X[apply(dd, 2, which.min), ]
      segments(x0=Ucap[,1], y0=Ucap[,2], x1=closest[,1], y1=closest[,2])
    }, silent = TRUE)
    if(inherits(tryPlot, "try-error"))
      tryPlotError(tryPlot)
  }

  return(list(
    # ............... arguments input ..........................
    X=X, Ntotal=Ntotal, sigma=sigma, beta1=beta1, nsurveys=nsurveys, xlim=xlim, ylim=ylim,
    # ............... generated values ..........................
    Habitat=as.vector(x), # a vector for the habitat covariate for each pixel
    Habgrid=gr,           # a 2-column matrix with the coordinates of each pixel
    nPix=nPix,            # the number of pixels in the study area
    N = N,                # true number of individuals per pixel
    U = U,                # locations of each individual in the population
    Ucap = Ucap,          # locations of each individual detected at least once
    nind=nind,            # the number of individuals detected at least once
    pixel=pixel)          # a matrix with a column for each survey and a row
        # for each individual detected at least once, with the pixel ID for the
        # activity centre or NA if the individual was not detected on the survey
  )
}

