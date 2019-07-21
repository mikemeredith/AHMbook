

# function by Mathias Tobler which computes the correlation matrix in residual occurrence in an latent-variable multi-species occupancy or N-mixture model as showcased on Chapter 8 of AHM2. Input must be the sims.list of the latent variables (LV), the sims.list of their coefficients (lv.coef) and the number of species. The function returns the residual correlation matrix, as described in Tobler et al. (Ecology, 2019) and other recent JSDM papers.

getLVcorrMat <- function(lv.coef, type=c("occupancy", "Nmix"), stat=mean){
  type <- match.arg(type)
  niter <- dim(lv.coef)[1]
  nspec <- dim(lv.coef)[2]
  cm.all <- array(NA, dim = c(niter, nspec, nspec))
  if(type == "occupancy") {
    eps.res <- apply(lv.coef, c(1,2), function(x) 1 - sum(x^2))
    fix <- diag(apply(eps.res, 2, mean))
  } else {
    fix <- 0
  }
  for(i in 1:niter){ # for each mcmc sample
    cm.all[i,,] <- cov2cor(tcrossprod(lv.coef[i,,]) + fix)
  }
  if(is.null(stat)) {
    return(cm.all)
  } else {
    return(apply(cm.all, c(2, 3), stat))
  }
}

# Function by Marc, email 2019-05-31

getEcorrMat <- function(beta, stat=mean){
  nspec <- dim(beta)[2]
  ecorraw <- apply(beta, 1, function(x) cor(t(x))) #  matrix, nspec^2 x niter
  ecorarr <- array(t(ecorraw), dim=c(dim(beta)[1], nspec, nspec))
  if(is.null(stat)) {
    return(ecorarr)
  } else {
    return(apply(ecorarr, c(2, 3), stat))
  }
}
