# Functions for the book Applied Hierarchical Modeling in Ecology (AHM)
# Marc Kery & Andy Royle, Academic Press, 2016.

# Function returning three fit-statistics (used in parboot GOF tests throughout book)
# (used, among others, in Chapter 7, e.g., Section 7.5.4)
fitstats <- function(fm) {
   observed <- unmarked::getY(fm@data)
   expected <- fitted(fm)
   resids <- residuals(fm)
   sse <- sum(resids^2, na.rm = T)                   # Sums of squares
   chisq <- sum((observed - expected)^2 / expected, na.rm = T) # Chisq
   freeTuke <- sum((sqrt(observed) - sqrt(expected))^2, na.rm = T) # F-T
   out <- c(SSE=sse, Chisq=chisq, freemanTukey=freeTuke)
   return(out)
}


# Define new fitstats function
#   (introduced in Section 7.9.3)
fitstats2 <- function(fm) {
   observed <- unmarked::getY(fm@data)
   expected <- fitted(fm)
   resids <- residuals(fm)
   n.obs <- apply(observed,1,sum,na.rm=TRUE)
   n.pred <- apply(expected,1,sum,na.rm=TRUE)
   sse <- sum(resids^2,na.rm=TRUE)
   chisq <- sum((observed - expected)^2 / expected,na.rm=TRUE)
   freeTuke <- sum((sqrt(observed) - sqrt(expected))^2,na.rm=TRUE)
   freeTuke.n<- sum((sqrt(n.obs)-sqrt(n.pred))^2,na.rm=TRUE)
   sse.n <- sum( (n.obs -n.pred)^2,na.rm=TRUE)
   chisq.n <- sum((n.obs - n.pred)^2 / expected,na.rm=TRUE)

   out <- c(SSE=sse, Chisq=chisq, freemanTukey=freeTuke,
      SSE.n = sse.n, Chisq.n = chisq.n, freemanTukey.n=freeTuke.n)
   return(out)
}

