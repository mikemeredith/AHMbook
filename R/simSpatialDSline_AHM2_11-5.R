
# Function for AHM2 Chapter 11, section 11.5

simSpatialDSline <-
function(N=1000, beta = 1, sigma=0.25, alpha0 = -2, W=1/2, L = 4, perp=FALSE, show.plots=TRUE){

  #   N: total population size in the rectangle
  #   beta: coefficient of SOMETHING on spatial covariate x
  #   sigma: scale of half-normal detection function
  #   W : truncation distance = strip half-width.

  # Create coordinates rasterized transect
  delta<- 0.1                 # '2D bin width'
  # Following code creates coordinates in order of "raster"
  #   (not necessary here but helps keep things organized)
  gry <- seq(delta/2, W*2 - delta/2, delta)
  ny <- length(gry)
  grx <- seq(delta/2, L - delta/2, delta)
  nx <- length(grx)
  grx <- rep(grx, ny)
  gry <- rev(sort(rep(gry,nx)))
  gr <- cbind(grx,gry)

  # Create spatially correlated covariate x and plot it
  V <- exp(-e2dist(gr,gr)/1)
  x <- t(chol(V))%*%rnorm(nrow(gr))
  r <- rasterFromXYZ(cbind(gr,x))

  # Simulate point locations as function of habitat covariate x
  probs <- exp(beta*x)/sum(exp(beta*x)) # probability of point in pixel (sum = 1)
  pixel.id <- sample(1:nrow(gr), N, replace=TRUE, prob=probs)
  # could simulate randomly within the pixel but it won't matter so place centrally
  u1 <- gr[pixel.id,1]
  u2 <- gr[pixel.id,2]
  # points(u1, u2, pch=20, col='black', cex = 1)  # plot points ### some pixels have >1 animal
  # title("Transect HDS")

  line.pts <- seq(0.01, L-0.01, .02)
  d.to.trap <- trap.ind <- rep(NA, N)
  traps <- cbind(line.pts, 0.5) # trap locations, spaced at 0.02 all along the line

  # Put together first part of output list
  outbasic <- list(
      # --------- arguments entered -------------------
      N=N, beta=beta, sigma=sigma, alpha0=alpha0, W=W, L=L,
      # --------- values generated ---------------------
      delta=delta,  # distance between pixel centres (spatial resolution of grid)
      grid=gr,      # 2-column matrix with x/y coordinates of all pixels
      Habitat=x,    # value of habitat covariate for each pixel
      Habraster = r,# a raster with the habitat covariate
      u1=u1, u2=u2, # x and y coordinates of all animals in the population
      traps=traps)  # 2-column matrix of trap locations


  if(!perp){
    dmat <- e2dist(cbind(u1,u2), traps)
    pbar <- trap.ind

    for(i in 1:nrow(dmat)){   # nrow(dmat) == N, ie, loop over animals
      ## This needs to loop over "traps" and flip a coin until encounter!
      haz <- exp(alpha0)*exp( -(dmat[i,]^2)/(2*sigma*sigma))
      probs <- 1-exp(-haz)
      captured <- rbinom(nrow(traps), 1, probs)
      pbar[i] <- 1- exp(-sum(haz))  # prob that animal will be detected at least once.
      if(sum(captured)==0)  # not captured, NAs everywhere
        next
      trap.ind[i] <- which(captured == 1)[1] # point on the line where animal first detected
      d.to.trap[i]<- dmat[i, trap.ind[i]]
    }
    # Trap.ind is where on the line guy was detected
    data <- cbind(trap.ind, d.to.trap, u1, u2)
    data <- data[!is.na(trap.ind),] # remove animals not detected
    pixel <- pixel.id[!is.na(trap.ind)] # pixels of animals detected

    outextra <- list(
      data=data,   # matrix with rows for each animal captured and columns for
              # trap of first capture, distance to trap, x and y coordinates.
      pbar=pbar,   # probability that the animal is captured at least once.
      pixel=pixel) # pixel ID of animals captured
  } else {
    # simulate ordinary DS data
    dmat <-  abs(u2 - 0.5)
    probs <- exp( -dmat*dmat/(2*sigma*sigma) )
    captured <- rbinom(N, 1, probs)

    data<- cbind(u1,u2)[captured==1,]
    outextra <- list(
      data=data, # a 2-column matrix with x and y coordinates of each animal captured.
      pixel = pixel.id[captured==1]) # pixel ID for each animal captured.
  }

  if(show.plots) {
    oldpar <- par(mar=c(3,3,3,6), "mfrow") ; on.exit(par(oldpar))
    tryPlot <- try( {
      image(r, col=topo.colors(10))
      abline(0.5, 0, lwd=2)
      image_scale(x, col=topo.colors(10))
      points(u1, u2, pch=20, col='black', cex = 1)  # plot points ### some pixels have >1 animal
      title("Transect HDS")

      if(perp) {
        for(i in 1:N){
          if(captured[i]==1) {
            points(u1[i], u2[i], col='red', cex=1.5)  # circle captured animals
            lines(c(u1[i], u1[i]), c(u2[i], 0.5) )    # lines for those captured
          }
        }
      } else {
          for(i in 1:nrow(dmat))
            lines(c(u1[i], traps[trap.ind[i],1]), c(u2[i], traps[trap.ind[i],2]) )
      }
    }, silent = TRUE)
    if(inherits(tryPlot, "try-error"))
      tryPlotError(tryPlot)
  }

  return(c(outbasic, outextra))
}
