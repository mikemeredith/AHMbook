# Functions for the book Applied Hierarchical Modeling in Ecology (AHM)
# Marc Kery & Andy Royle, Academic Press, 2016.

# wigglyOcc - AHM1 section 10.14 p622

# Function to generate a static occupancy data set with really wiggly covariate relationships
# in occupancy and detection probability
#   (introduced in AHM1 Section 10.14)
wigglyOcc <- function(seed = 1, show.plot = TRUE, verbose = TRUE){
# Function simulates really wiggly static site-occupancy data
#
# seed is for random number generator
#
# Choose sample sizes and seed for repeatability
M <- 240         # Number of sites
J <- 3           # Number of replicates
set.seed(seed)   # Allow repeatability

# Ecological process: Generation of latent occurrence state
Xsite <- seq(-2, 2, length.out = M)
psi <- c(seq(0.1, 0.9,,80), rep(0.9,,80), rep(0.3,,80))
z <- rbinom(M, 1, psi)

# Observation process: Generation of observed data
# Put covariate Xsurvey and p in order
Xsurvey <- seq(-2, 2,, M*J)
p.bp <- c(0, 0.6, 0.2, 0.9, 0.2, 0, 0.2) # "break points" for p model
p.ordered <- c(seq(p.bp[1], p.bp[2],, 120),
               seq(p.bp[2], p.bp[3],, 120),
               seq(p.bp[3], p.bp[4],, 120),
               seq(p.bp[4], p.bp[5],, 120),
               seq(p.bp[5], p.bp[6],, 120),
               seq(p.bp[6], p.bp[7],, 120))
x.index <- sample(1:length(Xsurvey))
Xsurvey <- matrix(Xsurvey[x.index], M, J, byrow = FALSE)
p <- matrix(p.ordered[x.index], M, J, byrow = FALSE)

# Sample detection/nondetection data
y <- array(dim = c(M, J))
for(j in 1:J){
   y[,j] <- rbinom(M, 1, z * p[,j])
}

# Look at data and produce some summaries
head(cbind(z = z, p = p, y = y))
if(verbose) {
  cat("   True number of occupied sites:", sum(z), "\n")
  cat("   Observed number of occupied sites:", sum(apply(y,1,max)), "\n")
  cat("   Proportional underestimation of distribution:", round((sum(z)-sum(apply(y,1,max)))/ sum(z), 2), "\n")
}

# Plot system (state and observation)
if(show.plot) {
  op <- par(mfrow = c(1,2), cex.main = 0.8) ; on.exit(par(op))
  tryPlot <- try( {
    plot(Xsite, psi,
        main = "Occupancy probability (red) and \nrealized presence/absence (black circles)",
        type = "l", ylim = c(-0.1, 1.1), col = "red", xlab = "Site covariate (Xsite)",
        lwd = 2, frame = FALSE)
    points(Xsite, jitter(z, amount = 0.02))
    plot(Xsurvey[order(x.index)], p[order(x.index)], type = "l", col = "red",
        main = "Detection probability (red) and \nobserved data (black circles)",
        xlab = "Survey covariate (Xsurvey)", ylab = "p", ylim = c(-0.1,1.1), lwd = 2, frame = FALSE)
    points(Xsurvey, jitter(y, amount = 0.02))
  }, silent = TRUE)
  if(inherits(tryPlot, "try-error"))
    tryPlotError(tryPlot)
}

return(list(M = M, J = J, Xsite = Xsite, Xsurvey = Xsurvey, psi = psi, z = z, p = p, y = y,
  x.index=x.index, p.ordered=p.ordered))
}

