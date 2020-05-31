
# Mizel, J.D., Schmidt, J.H., & Lindberg, M.S. (2018) Accommodating temporary emigration in spatial distance sampling models. Journal of Applied Ecology, 55, 1456-1464.
# Appendix S1. Simulation and JAGS code for the TPP model.
# Simulation code adapted from Kery and Royle (2016)
###############################################################################

# Changes to mesh with AHMbook conventions:
#   name sim.spatialHDS.TE -> sim.spatialDSte
#   b1 -> beta, int.lam -> lam0, adj.sigma ~~> sigma, T -> nsurveys
#   argument delta removed (values other than 1 throw an error)

simSpatialDSte <- function(
    nsites=28,     # number of sites
    dim=10,         # number of pixels along each side of the square site
    beta=1,         # the effect of habitat on the number of individuals in a pixel
    lam0=2.5,       # expected population size in the square
    nsurveys=4,     # number of surveys
    sigma=3,        # scale of half-normal detection function in pixels
    phi=0.6,        # availability
    theta=2,        # exponential correlation in the spatial covariate
    show.plots=3) {

  # Checks and fixes for input data -----------------------------
  nsites <- round(nsites[1])
  dim <- round(dim[1])
  stopifNegative(lam0, allowZero=FALSE)
  nsurveys <- round(nsurveys[1])
  stopifNegative(sigma, allowZero=FALSE)
  stopifnotProbability(phi)
  stopifNegative(theta, allowZero=FALSE)
  # --------------------------------------------


  npixels <- dim * dim
  B <- dim/2

  # sigma <- adj.sigma*B  # Default adj.sigma is 0.6 x radius B

  # Create coordinates for npixels x npixels grid
  delta <- 1
  grx <- seq(delta/2, 2*B - delta/2, delta) # mid-point coordinates
  gr <- expand.grid(grx,grx)         # Create grid coordinates
  center <- matrix(B,nrow=1,ncol=2)

  tr <- cbind(rep(B,length(grx)),grx)
  d1 <- e2dist(tr,gr)
  d <- apply(d1,2,min)

  # V <- exp(-e2dist(gr,gr)/1)
  # V <- exp(-e2dist(gr,gr)/2) #### changed 2019-05-20 v.0.1.4.9063
  V <- exp(-e2dist(gr,gr)/theta) #### changed 2019-05-20 v.0.1.4.9064

  # Create spatially correlated covariate x and plot it
  beta0 <- log(lam0/npixels) # intercept of log(N) ~ beta0 + beta(Habitat)
  x <- probs <- array(NA,dim=c(npixels,nsites))
  M <- rep(NA,nsites)

  for (j in 1:nsites){
    # z <-  t(chol(V))%*%rnorm(npixels)
    # x[,j]<- z
    x[,j] <- t(chol(V))%*%rnorm(npixels)     # habitat covariate for this site
    M[j] <- rpois(1, sum(exp(beta0 + beta*x[,j]))) # number of individuals at site ??
  }

  for (i in 1:npixels){
    for (j in 1:nsites){
      probs[i,j] <- exp(beta*x[i,j])/sum(exp(beta*x[,j])) # prob that animal at site j is at pixel i
    }
  }
  all.equal(colSums(probs), rep(1, nsites))

  Mind <- max(M)

  superpop <- array(0, c(Mind, nsites))
  for (j in 1:nsites){
    # ifelse(M[j]>0, superpop[1:M[j],j] <-1, superpop[,j] <- 0)
    superpop[1:M[j], j] <- M[j] > 0
  }
  all(colSums(superpop) == M)

  # Simulate individual locations for each survey
  #  All individuals get locations, ignore temp emigration at this stage
  pixel.id  <- array(NA, dim=c(Mind, nsites, nsurveys))
  for (i in 1:Mind){
    for (j in 1:nsites){
      for (k in 1:nsurveys){
        # pixel.id[i,j,k]  <-  sample(1:npixels, 1, replace=TRUE, prob=probs[,j])
        pixel.id[i,j,k]  <-  sample.int(npixels, 1, prob=probs[,j])
      }
    }
  }

  # y1 <- p <- array(NA,dim=c(Mind, nsites, nsurveys)) # p not used outside the loop
  y1 <- array(NA, dim=c(Mind, nsites, nsurveys))
  # Simulate observations
  # temp emigration recognised here -> availability
  # p = real member of superpop x availability x detection function (half-normal)
  for (i in 1:Mind){
    for (j in 1:nsites){
      for(k in 1:nsurveys){
        # p[i,j,k]<-superpop[i,j] * phi *
              # exp(-d[pixel.id[i,j,k]]*d[pixel.id[i,j,k]]/(2*(sigma^2)))
        # y1[i,j,k]<-rbinom(1, 1, p[i,j,k])
        p <- superpop[i,j] * phi *
              exp(-d[pixel.id[i,j,k]]*d[pixel.id[i,j,k]]/(2*(sigma^2)))
        y1[i,j,k]<-rbinom(1, 1, p)
      }
    }
  }
  Counts <- apply(y1,2:3, sum, na.rm=TRUE)

  pixel.id[y1==0] <- 0  # zap pixel.id for animals not detected or not real individual

  # Re-shape individual data structure into counts in site x pixel x visit array
  y <- array(NA, dim = c(nsites, npixels, nsurveys),
    dimnames = list(NULL, c(1:npixels)))
  for(i in 1:nsites){
    for (k in 1:nsurveys){
      # y[i,,k] <- table(factor(paste(pixel.id[,i,k], sep = ""), levels = c(1:npixels)))
      y[i,,k] <- tabulate(pixel.id[,i,k], nbins=npixels)
    }
  }

  # Do some plots
  if(show.plots > 0) {
    show.plots <- min(show.plots, nsites)
    oldpar <- par(mar=c(1,1,3,1), oma=c(2,0,2,0), "mfrow")
    oldAsk <- devAskNewPage(ask = dev.interactive(orNone = TRUE))
    on.exit({par(oldpar) ; devAskNewPage(oldAsk)})

    for(i in 1:show.plots) {
      if(nsurveys < 3)
        par(mfrow = c(1,2))
      if(nsurveys > 2)
        par(mfrow = c(2,2))
      img <- rasterFromXYZ(cbind(gr, x[,i]))
      for(j in 1:min(nsurveys, 4)) {
        raster::plot(img, col=rampBYR(255), axes=FALSE, box=FALSE)
        title(main=paste("survey", j), line=0.2)
        points(gr[pixel.id[, i, j], , drop=FALSE], pch=16)
        segments(dim/2, 0, dim/2, 10, lwd=3, col='black')  # The transect line
      }
      title(main=paste("Site", i, ": True population =", M[i]), cex.main=1.5, line=0, outer=TRUE)
      if(nsurveys > 4)
        mtext(paste(nsurveys - 4, "more surveys not shown"), side=1, outer=TRUE)
    }
  }

  return(list(  # ---------- arguments supplied -----------
    nsites=nsites, dim=dim, beta=beta, lam0=lam0, nsurveys=nsurveys,
    sigma=sigma, phi=phi, theta=theta,
    # ------------ values generated ---------------------------
    npixels=npixels,  # number of pixels in each site (= dim^2)
    B=B,                # distance from line to edge of square (= dim/2)
    M=M,                # true number of individuals at each site
    d=d,                # perpendicular distance of each pixel from the line
    Habitat=x,          # pixels x sites, value of habitat covariate for each pixel
    y=y,                # sites x pixels x surveys, number of animals detected
    Counts=Counts))     # sites x surveys, number of animals detected (summed over pixels)
}
