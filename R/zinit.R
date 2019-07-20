
# Generate starting values for survival analysis in JAGS or WinBUGS

# AHM2 section 3.2.3

zinit  <- function(CH){
  CH <- round(as.matrix(CH)) # could also be a data frame

  f <- suppressWarnings(apply(CH, 1, function(x) min(which(x!=0)))) # occasion of first capture
  zinit <- array(NA, dim = dim(CH))
  for(i in 1:nrow(CH)){
    if(f[i] >= ncol(CH)) # first captured on last occasion (or never!)
      next
    zinit[i,(f[i]+1):ncol(CH)] <- 1
  }
  return(zinit)
}





