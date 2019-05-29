

# function by Mathias Tobler which computes the correlation matrix in residual occurrence in an latent-variable multi-species occupancy or N-mixture model as showcased on Chapter 19 of the book. Input must be the sims.list of the latent variables (LV), the sims.list of their coefficients (lv.coef) and the number of species. The function returns the residual correlation matrix, as described in Tobler et al. (Ecology, 2019) and other recent JSDM papers.

getLVcorrMat <- function(LV, lv.coef, stat=mean){
  nspec <- dim(lv.coef)[2]
  cm.all.1 <- cm.all.2 <- array(NA, dim = c(dim(lv.coef)[1], nspec, nspec))
  eps.res <- apply(lv.coef, c(1,2), function(x) 1 - sum(x^2))
  for(i in 1:dim(lv.coef)[1]){ # for each mcmc sample
    # cm.all.1[i,,] <- cov2cor(tcrossprod(lv.coef[i,,]) + diag(dim(lv.coef[i,,])[1]))
    cm.all.2[i,,] <- cov2cor(tcrossprod(lv.coef[i,,]) + diag(apply(eps.res, 2, mean)))
  }
  cm.est <- apply(cm.all.2, c(2, 3), stat)
  return(cm.est)
}
