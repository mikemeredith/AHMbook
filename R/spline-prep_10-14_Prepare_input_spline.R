# Functions for the book Applied Hierarchical Modeling in Ecology (AHM)
# Marc Kery & Andy Royle, Academic Press, 2016.

# spline.prep - AHM1 section 10.14 p623

# Function to prepare input for BUGS model when fitting a spline for a covariate
#   (introduced in AHM1 Section 10.14)

spline.prep <- function(cov, nknot = NA){
# Function chooses knots and creates design matrices for fixed and
#   random-effects parts of a spline model for a chosen covariate
# Based on code by Crainiceanu et al. (2005) and Zuur et al. (2012)
# Allows you to choose number of knots or else uses it by the rule
#   given in Crainiceanu et al. (2005)
# Prepares fixed part of covariate as a quadratic polynomial

# Determine number and position of knots
# ifelse(is.na(nknot),
# n.knots <- max(5, min(round(length(unique(cov))/4), 35)),
# n.knots <- nknot)
if(is.na(nknot)) {
  n.knots <- max(5, min(round(length(unique(cov))/4), 35))
} else {
  n.knots <- nknot
}
prob.tmp <- seq(0,1, length = n.knots + 2)
prob <- prob.tmp[-c(1, length(prob.tmp))]
knots <- quantile(unique(cov), probs = prob)

# Create design matrices for fixed and random effects
X <- cbind(rep(1, length(cov)), cov, cov^2) # Fixed-eff DM
Z.tmp <- (abs(outer(cov, knots, "-")))^3
omega.all <- (abs(outer(knots, knots, "-")))^3
svd.omega.all <- svd(omega.all)
sqrt.omega.all <- t(svd.omega.all$v %*% (t(svd.omega.all$u) * sqrt(svd.omega.all$d)))
Z <- t(solve(sqrt.omega.all, t(Z.tmp)))  # Rand. eff. DM

# Output
return(list(cov = cov, knots = knots, X = X, Z = Z))
}


