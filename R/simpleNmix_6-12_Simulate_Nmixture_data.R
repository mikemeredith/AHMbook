# Functions for the book Applied Hierarchical Modeling in Ecology (AHM)
# Marc Kéry & Andy Royle, Academic Press, 2016.

# simpleNmix - section 6.12 p298

# Function to generate Nmix data under a time-for-space substitution design
#   (introduced in Section 6.12)
# Define function to simulate such data

simpleNmix <- function(nyear = 12, nrep = 4, beta0 = 2, beta1 = 0.1, alpha0 = 0.5, alpha1 = -0.1, alpha2 = 1){
# Simple function simulates data under binomial N-mixture model where you have
# a single site that is survyed over 'nyear' primary sampling periods
# ('seasons', 'years'), within which there are 'nrep' secondary samples each
# alpha0, alpha1 are the logit-linear coefficients of detection (p) on Time
#    and on a survey-specific covariate such as temperature (temp).
# beta0 and beta1 are the log-linear coefficients of expected abundance
#   (lambda) on Time.

Time <- 1:nyear
temp <- matrix(runif(nyear*nrep, -2, 2), ncol = nrep)
N <- rpois(n = nyear, lambda = exp(beta0 + beta1 * Time))
C <- array(NA, dim = c(nyear, nrep))
p <- plogis(alpha0 + alpha1*Time + alpha2*temp)
for(j in 1:nrep){
   C[,j] <- rbinom(n = nyear, size = N, prob =p[,j])
}
par(mfrow = c(3, 2))
curve(function(x) exp(beta0 + beta1 * x), 1, nyear, main = "Expected abundance (lambda) over time", frame = F, lwd = 2, ylab = "lambda", xlab = "Time")
plot(Time, N, main = "Realized abundance (N) over time", frame = F)
curve(function(x) plogis(alpha0 +alpha1 * x), 1, nyear, main = "p over time", frame = F, lwd = 2, xlab = "Time", ylab = "p (at averate temp)")
matplot(Time, C, main = "Counts (C) over time", frame = F)
curve(function(x) plogis(alpha0 + alpha2 * x), -2, 2, main = "p vs. Temperature", frame = F, lwd = 2, xlab = "Temperature", ylab = "p (at start of study)")
matplot(temp, C, main = "Counts (C) over time", frame = F)

return(list(nyear=nyear, nrep=nrep, beta0=beta0, beta1=beta1, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, N=N, C=C, Time=Time, temp = temp, p = p))
}

