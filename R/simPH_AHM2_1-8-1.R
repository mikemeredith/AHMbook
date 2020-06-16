
## AHM2 section 1.8.1

simPH <- function(
  # --- Sample sizes and design stuff ---
  npop = 18,               # Number of populations
  nyears = 17,              # Number of years (seasons)
  nreps = 10,               # Number of surveys per year (season)
  date.range = 1:150,      # Dates over which surveys may be conducted
  # --- Parameters of among-year dynamics ---
  initial.lambda = 300,    # Poisson mean of initial population size
  gamma.parms = c(0, 0.3), # mean and sd of lognormal interannual productivity
  # --- Parameters of within-year dynamics ---
  mu.range = c(50, 80),    # Range of date of peak flight period
                           #   (varies by site and year)
  sigma.range = c(10, 20), # Range of sigma of normal phenology curve
                           #   (varies by year only)
  # --- Parameters of observation process ---
  p.range = c(0.4, 0.6),   # Range of detection probabilities
                           #   (varies by site, year and visit)
  # --- Switch for plotting ---
  show.plot = TRUE)        # whether to browse plots or not
                           #   (should be set to FALSE when running sims)
  {  # -------------------- Start of function code -----------------

  # Function generates (insect) counts under a variant of a
  # 'phenomenological model' of Dennis et al. (JABES 2016).
  #
  # Interannual population model is exponential population growth,
  # with Poisson initial abundance governed by initial.lambda and
  # annually varying growth rate (or productivity parameter) gamma
  #
  # Within-year dynamics is described by a Gaussian curve with date of
  # mean flight period mu (site- and year-specific) and
  # length of flight period sigma (only year-specific).
  #
  # Counts are made subject to a detection probability (p), which varies
  # randomly according to a uniform distribution for every single count.
  #
  # Counts are plotted for up to 16 populations only.

  # Checks and fixes for input data -----------------------------
  npop <- round(npop[1])
  nyears <- round(nyears[1])
  nreps <- round(nreps[1])
  stopifnotInteger(date.range)
  stopifNegative(initial.lambda, allowZero=FALSE)
  stopifnotLength(gamma.parms, 2)
  stopifnotProbability(p.range)
  # ---------------------------------------------------------------

  # Simulate among-year population dynamics: exponential model for n
  n <- array(NA, dim = c(npop, nyears))  # Array for site-year abundance
  n[,1] <- rpois(npop, initial.lambda)
  gamma <- rlnorm(nyears-1, meanlog=gamma.parms[1], sdlog=gamma.parms[2])
  for(t in 2:nyears){
    n[,t] <- rpois(npop, n[,t-1] * gamma[t-1])
  }

  # Simulate within-year population dynamics: Normal curve for counts
  C <- date <- lambda <- a <- array(NA, dim = c(npop, nyears, nreps))  # Arrays for
  # site-year-visit counts, survey dates, relative pop. size and detection probability
  mu <- array(NA, dim = c(npop, nyears))  # Array for value of peak flight period date

  # Select survey dates, peak flight period (mu_it),
  # length of flight period sigma(t) and compute relative pop. size (a),
  # expected population size (lambda) and realized counts (C)

  # Draw annual value of flight period length (sigma)
  sigma <- runif(nyears, min(sigma.range), max(sigma.range))

  # Draw values of detection probability (p)
  p <- array(runif(prod(c(npop, nyears, nreps)), min(p.range), max(p.range)),
    dim = c(npop, nyears, nreps))

  # Compute and assemble stuff at the scale of the individual visit
  for(i in 1:npop){
    for(t in 1:nyears){
      # Survey dates for this yr and pop:
      survey.dates <- sort(round(
        runif(nreps, min(date.range), max(date.range))))
      date[i,t,] <- survey.dates     # Save these survey dates
      mu[i,t] <- runif(1, min(mu.range), max(mu.range)) # Flight peak
      for(k in 1:nreps){
        # a[i,t,k] <- (1 / (sigma[t] * sqrt(2 * pi)) ) * exp( -((date[i,t,k] - mu[i,t])^2) / (2 * sigma[t]^2) )             # Rel. population size
        a[i,t,k] <- dnorm(date[i,t,k], mu[i,t], sigma[t])
               # Rel. population size
        lambda[i,t,k] <- n[i,t] * a[i,t,k] * p[i,t,k] # Expected counts
        C[i,t,k] <- rpois(1, lambda[i,t,k])       # Realized counts
      }
    }
  }

  if(show.plot) {
    # Restore graphical settings on exit -------------------------
    oldpar <- par("mfrow", "mar")
    oldAsk <- devAskNewPage(ask = dev.interactive(orNone=TRUE))
    on.exit({par(oldpar); devAskNewPage(oldAsk)})
    # ------------------------------------------------------------

    # Simulate nice smooth normal curve for the plots
    nday <- length(date.range)
    aa <- ll <- array(NA, dim = c(npop, nyears, nday))  # Arrays
    pp <- array(runif(prod(c(npop, nyears, nday)), min(p.range), max(p.range)),
      dim = c(npop, nyears, nday))
    for(i in 1:npop){
      for(t in 1:nyears){
        for(k in 1:nday){
          aa[i,t,k] <- dnorm(date.range[k], mu[i,t], sigma[t])
            # Relative population size
          ll[i,t,k] <- n[i,t] * aa[i,t,k] * pp[i,t,k] # Expected counts
        }
      }
    }

    # Graphical output
    tryPlot <- try( {
      # Plot population dynamics and plot of all population sizes
      par(mfrow = c(2,1), mar = c(5,4,3,1))
      matplot(1:nyears, t(n), type = "l", lwd = 2, lty = 1,
          main = "Relative population size (n) for each population and year",
          ylab = "n", xlab = "Year", frame = FALSE, xaxt='n')
      tmp <- pretty(1:nyears)
      tmp[1] <- 1
      axis(1, at=tmp)
      plot(table(n), xlab = 'Relative population size', ylab = 'Frequency', 
          main = 'Frequency distribution of relative population size\nfor all sites and years',
          frame = FALSE)

      # Plot time-series of relative expected abundance for up to 16 populations
      par(mfrow = c(4,4), mar = c(5,4,3,1))
      limit <- ifelse(npop < 17, npop, 16)
      for(i in 1:limit){     # Plot only for 4x4 populations
        matplot(date.range, t(aa[i,,]), type = "l", lty = 1, lwd = 2,
            ylim = c(0, max(aa[i,,])), xlab = "Date", ylab = "Rel. abundance", 
            main = paste("Phenology in pop ", i, sep = ''), frame = FALSE)
      }

      # Plot time-series of relative expected abundance for up to 16 populations
      par(mfrow = c(4,4), mar = c(5,4,3,1))
      limit <- ifelse(npop<17, npop, 16)
      for(i in 1:limit){     # Plot only for 4x4 populations
        matplot(date.range, t(ll[i,,]), type = "l", lty = 1, lwd = 2,
            ylim = c(0, max(ll[i,,])), xlab = "Date", ylab = "Exp. abundance",
            main = paste("Rel. exp. n in pop ", i, sep = ''), frame = FALSE)
      }

      # Plot time-series of counts (= relative, realized abundance) for up to 16 populations
      par(mfrow = c(4,4), mar = c(5,4,3,1))
      limit <- ifelse(npop<17, npop, 16)
      for(i in 1:limit){     # Plot only for 4x4 populations
        matplot(t(date[i,,]), t(C[i,,]), type = "b", lty = 1, lwd = 2,
            ylim = c(0, max(C[i,,])), xlab = "Date", ylab = "Counts",
            main = paste("Pop ", i, "(mean n =", round(mean(n[i,])), ")"), frame = FALSE)
      }
    }, silent = TRUE)
    if(inherits(tryPlot, "try-error"))
      tryPlotError(tryPlot)
  }
  # Numerical output
  return(list(
    # ---------- arguments input --------------------------
    npop = npop, nyears = nyears, nreps = nreps, date.range = date.range,
    initial.lambda = initial.lambda, gamma.parms = gamma.parms,
    mu.range = mu.range, sigma.range = sigma.range, p.range = p.range,
    # ------------ generated values -----------------------
    # abundance
    gamma = gamma,   # nyears-1 vector, change in abundance
    n = n,           # site x year matrix, true abundance
    # phenology
    mu = mu,         # site x year matrix, mean of the flight period
    sigma = sigma,   # nyears vector, half-length of flight period
    # detection
    date = date,     # site x year x nreps, dates of the surveys
    a = a,           # site x year x nreps, phenology term
    lambda = lambda, # site x year x nreps, expected counts
    p = p,           # site x year x nreps, probability of detection
    C = C))          # site x year x nreps, simulated counts
}
# -------------------- End of function definition -----------------

