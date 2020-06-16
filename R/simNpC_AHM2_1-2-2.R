
# AHM2 12.2

# It generates counts from a single population observed over T years and which can be observed with or without imperfect detection, ... The goal of this function is to focus on what happens with relative-abundance inference when temporal patterns in abundance are confounded with temporal patterns in detection probability. Hence, we can simulate a stable population or one with linear increase or decrease with specified start and end points, and around which there is Poisson noise. The observed counts are Binomial outcomes with a detection probability which can similarly be chosen to be constant or change linearly over time.

# Define the simulation function
simNpC <- function(
  T = 20,            # length of time series
  expN = c(100, 75), # expected abundance at start and end of period, linear trend
  dp = c(0.5, 0.5),  # detection probability at start and end of period, linear trend
  show.plot = TRUE)  # whether to show plots or not
  {

  # Function simulates a single time-series of counts of length T
  # years from a population with a linear trend in expected abundance
  # (expN) leading from expN[1] to expN[2] and with a linear trend in
  # detection probability (dp) leading from dp[1] to dp[2].
  # Default is for a strongly declining population with constant p = 0.5.

  # ---- Checks and fixes for input data  -------------------
  T <- round(T[1])
  stopifnotLength(expN, 2)
  stopifNegative(expN)
  stopifnotLength(dp, 2)
  stopifnotProbability(dp)
  # -----------------------------------------------------
  
  # Pick values of expected abundance (lam) and
  #   detection probability (dp) for each year
  lambda <- seq(expN[1], expN[2], length.out = T)
  p <- seq(dp[1], dp[2], length.out = T)

  # Draw realized abundance (N) and the observed counts (C)
  N <- rpois(T, lambda)
  C <- rbinom(T, N, p)
  # Plots
  if(show.plot) {
    oldpar <- par(mfrow = c(1, 3), mar = c(5,5,1,1), cex.axis = 1.2,
      cex.lab = 1.2, cex = 1.2)
        on.exit(par(oldpar))
    tryPlot <- try( {
      plot(1:T, lambda, xlab = 'Year', ylab = 'Expected abundance (lambda)',
          ylim = c(0, max(expN)), type = 'l', lwd = 3, col = 2, frame = FALSE)
      plot(1:T, p, xlab = 'Year', ylab = 'Detection prob. (p)',
          ylim = c(0, 1), type = 'l', lwd = 3, col = 4, frame = FALSE)
      plot(1:T, N, xlab = 'Year', ylab = 'Counts, Abundance',
          ylim = c(0, max(N)), pch = 16, frame = FALSE)
      points(1:T, C, pch = 1)
      lines(1:T, lambda, col = 2, lwd = 2)
      lines(1:T, lambda*p, col = 1, lwd = 2, lty=2)
      legend(1, 0.24*max(N), c('True N', 'Observed C'), pch = c(16,1),
          cex = 0.8, bty = 'n')
      legend(1, 0.14*max(N), c('Expected N (lambda)',
          'Exp. relative abundance\n (lambda * p)'), lty = c(1,2), lwd = 3,
          col = c(2,1), cex = 0.8, bty = 'n')
    }, silent = TRUE)
    if(inherits(tryPlot, "try-error"))
      tryPlotError(tryPlot)
  }

  # Output
  return(list(
    # ---------- arguments input --------------------------
    T = T, expN = expN, dp = dp,
    # ------------ generated values -----------------------
    lambda = lambda,  # expected abundance for each year
    p = p,            # detection probability (dp) for each year
    N = N,            # realised abundance
    C = C))           # observed counts
}
