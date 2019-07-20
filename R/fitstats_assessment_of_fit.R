# Functions for the book Applied Hierarchical Modeling in Ecology (AHM)
# Marc Kery & Andy Royle, Academic Press, 2016.

# Function returning three fit-statistics (used in parboot GOF tests throughout book)
# (used, among others, in AHM1 Chapter 7, e.g., AHM1 Section 7.5.4)

# Updated 2019-01-14 to cope with NAs in the data, see AHM2 2.3.3

fitstats <- function(fm) {
  observed <- unmarked::getY(fm@data)
  notna <- !is.na(observed) # to accommodate missing values
  expected <- fitted(fm)
  resids <- residuals(fm)
  sse <- sum(resids[notna]^2)
  chisq <- sum((observed[notna] - expected[notna])^2 / expected[notna])
  freeTuke <- sum((sqrt(observed) - sqrt(expected))[notna]^2)
  out <- c(SSE=sse, Chisq=chisq, freemanTukey=freeTuke)
  return(out)
}


# Define new fitstats function
#   (introduced in AHM1 Section 7.9.3)
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

