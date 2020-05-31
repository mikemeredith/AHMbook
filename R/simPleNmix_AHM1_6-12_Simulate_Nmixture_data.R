# Functions for the book Applied Hierarchical Modeling in Ecology (AHM)
# Marc Kery & Andy Royle, Academic Press, 2016.

# simpleNmix - AHM1 section 6.12 p298

# Function to generate Nmix data under a time-for-space substitution design
#   (introduced in AHM1 Section 6.12)
# Define function to simulate such data

simpleNmix <- function(nyears = 12, nreps = 4, beta0 = 2, beta1 = 0.1, alpha0 = 0.5, alpha1 = -0.1, alpha2 = 1, show.plot = TRUE){
# Simple function simulates data under binomial N-mixture model where you have
# a single site that is survyed over 'nyears' primary sampling periods
# ('seasons', 'years'), within which there are 'nreps' secondary samples each
# alpha0, alpha1 are the logit-linear coefficients of detection (p) on Time
#    and on a survey-specific covariate such as temperature (temp).
# beta0 and beta1 are the log-linear coefficients of expected abundance
#   (lambda) on Time.
if(FALSE) x <- NULL  # Fix issues with 'curve'

# Checks and fixes for input data -----------------------------
nyears <- round(nyears[1])
nreps <- round(nreps[1])
# --------------------------------------------

Time <- 1:nyears
temp <- matrix(runif(nyears*nreps, -2, 2), ncol = nreps)
N <- rpois(n = nyears, lambda = exp(beta0 + beta1 * Time))
C <- array(NA, dim = c(nyears, nreps))
p <- plogis(alpha0 + alpha1*Time + alpha2*temp)
for(j in 1:nreps){
   C[,j] <- rbinom(n = nyears, size = N, prob =p[,j])
}

if(show.plot) {
  op <- par(mfrow = c(3, 2)) ; on.exit(par(op))
  curve(exp(beta0 + beta1 * x), 1, nyears, main = "Expected abundance (lambda) over time", frame = FALSE, lwd = 2, ylab = "lambda", xlab = "Time")
  plot(Time, N, main = "Realized abundance (N) over time", frame = FALSE)
  curve(plogis(alpha0 +alpha1 * x), 1, nyears, main = "p over time", frame = FALSE, lwd = 2, xlab = "Time", ylab = "p (at averate temp)")
  matplot(Time, C, main = "Counts (C) over time", frame = FALSE)
  curve(plogis(alpha0 + alpha2 * x), -2, 2, main = "p vs. Temperature", frame = FALSE, lwd = 2, xlab = "Temperature", ylab = "p (at start of study)")
  matplot(temp, C, main = "Counts (C) over time", frame = FALSE)
}
return(list(nyears=nyears, nreps=nreps, beta0=beta0, beta1=beta1, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, N=N, C=C, Time=Time, temp = temp, p = p))
}

