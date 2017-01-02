
# Function to extract the biggest cross-correlations from an mcmc or mcmc.list object

bigCrossCorr <- function(x, big = 0.6, digits=3) {
  # Function to extract the biggest cross-correlations
  #  from an mcmc or mcmc.list object.
  #
  # x : an mcmc or mcmc.list object as returned by rjags::jags or jagsUI::jags.samples
  # big : only values outside the range -big to +big will be returned
  # digits : number of digits to return.
  #
  # Returns a data frame with 3 columns, for the names of the two parameters and
  #  the correlation coefficient.
  # See ?coda::crosscorr for details
  #
  # Mike Meredith, 1 Jan 2017

  if(!inherits(x, c("mcmc", "mcmc.list")))
    stop("'x' must be an 'mcmc' or 'mcmc.list' object.")
  xcor <- coda::crosscorr(x)
  xcor[lower.tri(xcor, diag=TRUE)] <- 0
  BIG <- which(abs(xcor) > big, arr.ind=TRUE)
  nms <- rownames(xcor)
  return(data.frame(par1=nms[BIG[, 1]], par2=nms[BIG[, 2]],
    corr=round(diag(xcor[BIG[, 1], BIG[, 2]]), digits)))
}

