
# AHM2 Chapter 10

simPPe <- function(lscape.size = 150, buffer.width = 25, variance.X = 1, theta.X = 10, M = 250, beta = 1, quads.along.side = 6, show.plots = TRUE) {
  #
  # Name means 'SIMulate Point Pattern Educational version'

  # Function to simulate spatial point pattern
  # in a heterogeneous landscape simulated on a square grid.
  # The study area ('core') is simulated inside of a larger landscape that
  # includes a buffer. The size of the core is defined by the lscape.size
  # minus twice the buffer.
  # There is one habitat covariate X that affects the intensity
  # of the points.
  # X is spatially structured with negative exp. spatial autocorrelation;
  # the parameters of the field can be chosen to create large
  # 'islands' of similar values or no 'islands' at all,
  # in which case the field is spatially unstructured.
  #
  # The intensity of STATIC points (e.g. home-range centers)
  # may be inhomogeneous and affected by the coefficient beta,
  # which is the log-linear effect of X.
  #
  # *** Function arguments ***
  # lscape.size:  total side length of simulated landscape (length in units,
  #               this includes a buffer)
  # buffer.width: width of buffer around core study area
  # variance.X:   variance of Gaussian random field (covariate X)
  # theta.X:      scale parameter of correlation in field (must be >0)
  # M:            Expected number of activity centers in core area
  # beta:         coefficient of the habitat
  # quads.along.side: number of quadrats along the side of the core area
  #    (this is the parameter of the gridding process
  #    and determines the size of the quadrats)

  # -------------- Check and fix input -----------------------
  buffer.width <- round(buffer.width[1])
  stopifNegative(buffer.width)
  quads.along.side <- round(quads.along.side[1])
  stopifNegative(quads.along.side, allowZero=FALSE)
  lscape.size <- round(lscape.size[1])
  stopifnotGreaterthan(lscape.size, 2 * buffer.width + quads.along.side - 1)
  variance.X <- variance.X[1]
  stopifNegative(variance.X)
  theta.X <- theta.X[1]
  stopifNegative(theta.X)
  M <- round(M[1])
  stopifNegative(M)
  # ------------------------------------------------------------

  # --------------- Define basic geometry of simulation --------------
  #
  # Size of core study area (the 'core') and its proportion of total landscape area
  size.core <- lscape.size - 2 * buffer.width
  prop.core <- size.core^2/lscape.size^2 # ratio core area / total area

  # Discrete approximation of total landscape
  pixel.size <- 1      # length of side of square pixel of simulated landscape

  # Coordinates (mid-points of basic pixel unit of simulation) and lscape
  x <- seq(1, lscape.size, pixel.size)-0.5    # x coordinate of pixels
  y <- seq(1, lscape.size, pixel.size)-0.5    # y coordinate of pixels
  grid <- as.matrix(expand.grid(x,y))       # resulting grid

  # Compute lambda of point pattern: limit of expected number of points per areal unit, when latter goes towards zero
  lambda_pp <- M / size.core^2

  # Define a core area in the middle of the square
  # This core is then divided up to define a number of quadrats
  # within which abundance and occurrence are measured
  quad.size <- size.core / quads.along.side
  breaks <- seq(buffer.width, size.core+buffer.width, by = quad.size) # boundaries of quadrats
  mid.pt <- breaks[-length(breaks)] + 0.5 * quad.size # quadrat mid-points
  core <- range(breaks)   # range of x and y coordinates in the core
  nsite <- length(mid.pt)^2

  # Simulate habitat covariate: a spatially correlated Gaussian random field (i.e., a Gaussian random field with negative exponential corr.)
  if(requireNamespace("RandomFields", quietly=TRUE)) {
    RandomFields::RFoptions(seed=NA)
    field <- matrix(RandomFields::RFsimulate(RandomFields::RMexp(var = variance.X, scale = theta.X),
      x=x, y=y, grid=TRUE)@data$variable1, ncol = lscape.size) # MVN r.v. with spatial correlation
  } else {
    message("Using package 'fields' instead of 'RandomFields'; see help(simPPe).")
    obj <- circulantEmbeddingSetup(grid=list(x=x, y=y), Covariance="Exponential", aRange=theta.X)
    tmp <- try(circulantEmbedding(obj), silent=TRUE)
    if(inherits(tmp, "try-error"))
      stop("Simulation of random field failed.\nTry with smaller values for 'lscape.size' or 'theta.X'.")
    field <- matrix(tmp * sqrt(variance.X), ncol = lscape.size)
  }

  # --------------- Simulate points in the field --------------------
  #
  # Simulate binomial point process for activity centers
  M2 <- round(M/prop.core) # Number of individuals in the total landscape
                           #   (incl. buffer)

  # Simulate point locations as function of habitat covariate x
  probtemp <- exp(beta[1]*c(field)) # log-linear model for intensity on X
  probs <- probtemp / sum(probtemp) # normalize to get probability of getting a point in a pixel
  pixel.id <- sort(sample(1:lscape.size^2, M2 , replace=TRUE, prob=probs))

  # Simulate locations randomly within the pixel (unlike sim.spatialDS)
  u1 <- grid[pixel.id,1] + runif(M2, -pixel.size/2, pixel.size /2)
  u2 <- grid[pixel.id,2] + runif(M2, -pixel.size/2, pixel.size /2)
  u <- cbind(u1, u2)    # collect AC coordinates together


  # ------ Summarization of point pattern within quadrats -------------
  #
  # This is INSIDE of the observation core
  #
  # Summarization for abundance (N) at every quadrat
  Nac <- as.matrix(table(cut(u[,1], breaks=breaks),
     cut(u[,2], breaks= breaks))) # quadrat-specific abundance for AC
  E_N <- round(mean(Nac),2)       # average realized abundance per quadrat

  # Summarization for presence/absence (z) at every quadrat
  zac <- Nac   ;  zac[zac>1] <- 1 # quadrat-specific occurrence for AC
  E_z <- round(mean(zac), 2)      # proportion occupied quadrats

  # ------------------ Visualizations ---------------------------
  if(show.plots) {
    oldpar <- par(mfrow = c(1, 3), mar = c(4,2,5,2), cex.main = 1.8, cex.axis = 1.2) ; on.exit(par(oldpar))

    tryPlot <- try( {
      # *** Fig. 1: Original point pattern
      # Random field of X with activity-centers overlaid
      image(rasterFromXYZ(cbind(grid, c(field))), col=topo.colors(10),
          main = "Point pattern with\ncore and buffer area",
          xlab = "", ylab = "", axes = FALSE, asp = 1)
      mtext(paste("Mean intensity (lambda) =", round(lambda_pp, 5)), side=1)
      polygon(c(buffer.width, size.core+buffer.width, size.core+buffer.width, buffer.width),
            c(buffer.width, buffer.width, size.core+buffer.width, size.core+buffer.width), 
            lwd = 2, lty = 1)
      points(u[,1], u[,2], pch=20, col='black', cex = 1.2)  # plot points
      # points(u1, u2, pch=20, col='black', cex = 1.2)  # plot points

      # *** Fig. 2: Show abundance and presence/absence in each quadrat on original landscape ***
      # Covariate 1: the Gaussian random field with autocorrelation
      # Reproduce random field with activity centers
      image(rasterFromXYZ(cbind(grid, c(field))), col=topo.colors(10), main = "Abundance, N",
          xlab = "", ylab = "", axes = FALSE, asp = 1)
      mtext(paste0("Mean(N) = ", E_N, ", var(N) = ", round(var(c(Nac)), 2)), side=1)
      polygon(c(buffer.width, size.core+buffer.width, size.core+buffer.width, buffer.width),
          c(buffer.width, buffer.width, size.core+buffer.width, size.core+buffer.width), lwd = 2, lty = 1)
      # Add activity centers
      points(u[,1], u[,2], pch=20, col='black', cex = 1.2)  # plot points
      # Overlay survey quadrats
      for(i in 1:length(breaks)){
         for(k in 1:length(breaks)){
           segments(breaks[i], breaks[k], rev(breaks)[i], breaks[k])
           segments(breaks[i], breaks[k], breaks[i], rev(breaks)[k])
         }
      }
      # Print abundance into each quadrat
      for(i in 1:length(mid.pt)){
        for(k in 1:length(mid.pt)){
          text(mid.pt[i], mid.pt[k], Nac[i,k], cex =4^(0.8-0.5*log10(quads.along.side)), col='red')
        }
      }

      # Figure 3 for presence/absence of activity centers (= distribution)
      # Reproduce random field with activity centers
      image(rasterFromXYZ(cbind(grid, c(field))), col=topo.colors(10), main = "Occurrence, z",
          xlab = "", ylab = "", axes = FALSE, asp = 1)
      mtext(paste("Mean(z) =", E_z), side=1)
      polygon(c(buffer.width, size.core+buffer.width, size.core+buffer.width, buffer.width),
          c(buffer.width, buffer.width, size.core+buffer.width, size.core+buffer.width), lwd = 2, lty = 1)
      # Add activity centers
      points(u[,1], u[,2], pch=20, col='black', cex = 1.2)  # plot points
      # Overlay quadrats
      for(i in 1:length(breaks)){
         for(k in 1:length(breaks)){
           segments(breaks[i], breaks[k], rev(breaks)[i], breaks[k])
           segments(breaks[i], breaks[k], breaks[i], rev(breaks)[k])
         }
      }
      # Print presence/absence into each quadrat
      for(i in 1:length(mid.pt)){
        for(k in 1:length(mid.pt)){
          text(mid.pt[i], mid.pt[k], zac[i,k], cex =4^(0.8-0.5*log10(quads.along.side)), col='red')
        }
      }

      # Mike: Shade UNoccupied quadrats (which have abundance N = 0 or occurrence z = 0)
      for(i in 1:(length(breaks)-1)){
        for(k in 1:(length(breaks)-1)){
          if(zac[i,k] == 1) # grey-out UNoccupied quads
            next
          polygon(c(breaks[i], breaks[i+1], breaks[i+1], breaks[i]),
              c(breaks[k], breaks[k], breaks[k+1], breaks[k+1]),
              col = adjustcolor("black", 0.6))
        }
      }
    }, silent = TRUE)
    if(inherits(tryPlot, "try-error"))
      tryPlotError(tryPlot)
  }
  
# Numerical output
return(list(
  # ----------------- arguments input -----------------------
  grid.size = lscape.size, buffer.width = buffer.width, variance.X = variance.X, theta.X = theta.X, M = M, beta = beta, quads.along.side = quads.along.side,
  # ---------------- generated values -------------------------
  core = core,            # range of x and y coordinates in the 'core'
  M2 = M2,                # Number of ACs in the total landscape (incl. buffer)
  grid = grid,            # Coordinates of the centre of each pixel.
  pixel.size = pixel.size,# 1; length of side of square pixel of landscape
  size.core = size.core,  # the width = height of the core area
  prop.core = prop.core,  # proportion of the landscape in the core
  X = field,              # lscape.size x lscape.size matrix of covariate values for each pixel
  probs = probs,          # corresponding matrix of probability of AC in pixel (sums to 1)
  pixel.id = pixel.id,    # M2 vector, which pixel each AC is inside.
  u = u,                  # M2 x 2 matrix, coordinates of each AC
  nsite = nsite,          # total number of quadrats in the core
  quad.size = quad.size,  # width = height of each quadrat
  breaks = breaks,        # boundaries of each quadrat
  mid.pt = mid.pt,        # mid=points of each quadrat
  lambda_pp = lambda_pp,  # intensity of point pattern (ACs per unit area)
  Nac = Nac,    # matrix, quads.along.side x quads.along.side, site-specific abundance of ACs
  zac = zac,    # matrix, quads.along.side x quads.along.side, 0/1 occurrence
  E_N = E_N,    # scalar, average realized abundance per quadrat.
  E_z = E_z))   # scalar, average realized occupancy per quadrat.
} # ------------ End of function definition --------------------
