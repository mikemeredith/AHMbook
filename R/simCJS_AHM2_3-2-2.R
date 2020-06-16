# AHM2 section 3.2.2

# ------ Start of definition of the data simulation function ------
simCJS <- function(
    n.occ = 6,         # number of occasions (e.g., years)
    n.marked = 20,     # number of marked individuals per occasion
    phi = 0.7,         # apparent survival probability
    p = 0.4,           # recapture probability
    show.plot = TRUE)  # whether to show plots or not
  {  # -------- start of function code ----------------------
  # This function generates individual capture-histories under a
  # CJS model with possibly time-dependent parameters. It is based
  # on code written by Michael Schaub for chap. 7 of the BPA book.
  # The number of values for interval-specific survival (phi) and
  # time-specific detection (p) must be ensured to be equal to the
  # number of occasions (n.occ) minus 1
  #
  # Written by Marc Kery, Sep 2018
  #
  # Changes by Mike Meredith:
  #  Allow n.marked, phi and p to be EITHER scalar OR vector of length n.occ-1
  #  ... and catch errors.
  #  Modified code for generating z and ch (gives different values with set.seed)
  #  Restore ask and par on exit.

  # Check and fixes for input data -------------------------
  n.occ <- round(n.occ[1])
  n.marked <- round(n.marked)
  stopifnotLength(n.marked, n.occ-1, allow1=TRUE)
  stopifnotLength(phi, n.occ-1, allow1=TRUE)
  stopifnotProbability(phi)
  stopifnotLength(p, n.occ-1, allow1=TRUE)
  stopifnotProbability(p)
  # --------------------------------------------------------

  # Deal with input
  if(length(n.marked) == 1)
    n.marked <- rep(n.marked, n.occ-1)   # Annual number of newly marked individuals
  n.ind <- sum(n.marked)
  if(length(phi) == 1)
    phi <- rep(phi, n.occ-1)
  if(length(p) == 1)
    p <- rep(p, n.occ-1)

  # Vector (f) with marking occasion (ie, first capture occasion)
  f <- rep(1:length(n.marked), n.marked)

  # Fill the true state matrix (z) and capture-history matrix (ch)
  z <- matrix(NA, nrow = n.ind, ncol = n.occ) # true states z
  ## Mike says: easier to remove NAs in ch than to insert them into z
  ch <- matrix(0, nrow = n.ind, ncol = n.occ) # observed capture history
  z[f == 1, 1] <- 1   # animals caught on first occasion definitely alive
  ch[f == 1, 1] <- 1  # ... and definitely caught
  for(t in 2:n.occ)  {
    z[, t] <- suppressWarnings(rbinom(n.ind, 1, z[, t-1]*phi[t-1])) # have they survived?
    ## Mike says: rbinom gives NA if z is NA (and warns); that's what we want here
    ch[, t] <- suppressWarnings(rbinom(n.ind, 1, z[, t]*p[t-1]))    # were they caught?
    ## Mike says: NA is not what we want here, but easy to fix later
    z[f == t, t] <- 1   # animals first caught on occasion t definitely alive
    ch[f == t, t] <- 1  # ... and definitely caught
  }
  ch[is.na(ch)] <- 0 # fix the unwanted NAs

  # Tally up number alive, marked and in study area
  n.alive <- colSums(z, na.rm = TRUE) ## Mike is a fan of colSums and rowSums!

  # Visualizations
  if(show.plot){
    # Restore graphical settings on exit
    oldpar <- par(mfrow = c(1, 1), mar = c(5,5,5,3), cex.lab = 1.3, cex.axis = 1.3)
    oldAsk <- devAskNewPage(ask = dev.interactive(orNone=TRUE))
    on.exit({par(oldpar); devAskNewPage(oldAsk)})

    tryPlot <- try( {
      # PLOT 1
      # Plot trajectory of phi and p
      plot(1:(n.occ-1), phi, typ= 'n', ylim = c(0, 1),
        frame = FALSE, main = 'Trajectories of phi and p')
      points(1:(n.occ-1), phi, type= 'b', cex = 2, pch = 16, col = 2)
      points(1:(n.occ-1), p, type= 'b', cex = 2, pch = 16, col = 4, lty=3)
      legend('top', legend = c('Apparent survival (phi)', 'Recapture (p)'),
        lty=c(1,3), lwd=2, col=c(2,4), pch=16, cex=2,
        inset=c(0, -0.05), bty='n', xpd=NA, horiz=TRUE)

      # PLOT 2
      par(mfrow = c(2, 2))
      # Plot the true alive/dead pattern (z)
      mapPalette <- colorRampPalette(c("white", "black"))
      image(x = 1:n.occ, y = 1:n.ind, z = t(z), col = mapPalette(10), axes = TRUE,
        xlab = "Year", ylab = "Individual",
        main = 'z matrix of latent states in the CJS model: \nAlive (black) or dead (white) per individual and occasion')

      # Plot the observed alive/dead pattern (y, or ch)
      image(x = 1:n.occ, y = 1:n.ind, z = t(ch), col = mapPalette(10), axes = TRUE,
        xlab = "Year", ylab = "Individual",
        main = 'Observed data = capture-history matrix ch in the CJS model: \nDetected (black) or not detected (white) per individual and occasion')
      box()

      # Superimpose the two images
      tmp <- z    # copy z into tmp
      tmp[z==1 & ch == 0] <- -1.1  # Mark detection errors as -1
      mapPalette <- colorRampPalette(c("blue", "white", "black"))
      image(x = 1:n.occ, y = 1:n.ind, z = t(tmp), col = mapPalette(10), axes = TRUE,
      xlab = "Year", ylab = "Individual", main = 'Combopic of z and ch: not in study (white), alive & detected (black), \nalive & undetected (blue) and dead (grey) per individual and occasion')

      # Population size trajectory of marked and alive in study area
      plot(1:n.occ, n.alive, xlab = 'Year', ylab = 'Number alive',
      main = 'Number of marked animals alive and in study area', frame = FALSE,
      type = 'b', cex = 2, pch = 16, ylim = c(0, max(n.alive)))
    }, silent = TRUE)
    if(inherits(tryPlot, "try-error"))
      tryPlotError(tryPlot)
  }

  # Output
  return(list(
    # ---------- arguments input --------------------------
    n.occ = n.occ, n.marked = n.marked, phi = phi, p = p,
    # ------------ generated values -----------------------
    z = z,              # n.ind x n.occ matrix, 1 if alive and in study area
    ch = ch,            # n.ind x n.occ matrix, 1 if captured
    f = f,              # n.ind vector, occasion marked (= first capture)
    n.ind = n.ind,      # scalar, total number of individuals marked
    n.alive = n.alive)) # n.occ vector, number alive and in study area
}
# ------ End of definition of the data simulation function ------

