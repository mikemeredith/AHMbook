

# Function to estimate total population size of Great tits
#   in Switzerland from unmarked fit (introduced in Section 6.9.4)
Nhat <- function(fm = fm5ZIP, iLength = 0, area = 1) {
   betavec <- coef(fm)[1:8]
   DM <- cbind(rep(1,length(pelev)), pelev, pelev^2, pforest,
   pforest^2, rep(iLength,length(pelev)), pelev*pforest, pelev*pforest^2)
   pred <- exp(DM %*% betavec) * (1-plogis(coef(fm)[25])) * (1/area)
   pred2 <- pred
   N1 <- sum(pred[-which(CH$elev > 2250),])# Drop quads > 2250 m
   pred2[pred2 > 100] <- 100      # Censor freak-high estimates
   N2 <- sum(pred2[-which(CH$elev > 2250),]) # Estimate with freaks censored
   out <- which(CH$water > 50 | CH$elev > 2250)
   N3 <- sum(pred2[-out,])  # Estimate excluding water bodies
   return(c(N1 = N1, N2 = N2, N3 = N3))
}

