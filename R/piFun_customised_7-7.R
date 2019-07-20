# Functions for the book Applied Hierarchical Modeling in Ecology (AHM)
# Marc Kery & Andy Royle, Academic Press, 2016.

# instRemPiFun, crPiFun, crPiFun.Mb, MhPiFun  - AHM1 section 7.7 p346
# Customised 'piFun's for unmarked::multinomPois

# pi function for removal design with 3 intervals of unequal length (2, 3, 5 minutes)
#    (introduced in AHM1 Section 7.7)
instRemPiFun <- function(p){
   M <- nrow(p)
   J <- ncol(p)
   pi <- matrix(NA, M, J)
   p[,1] <- pi[,1] <- 1 - (1 - p[,1])^2
   p[,2] <- 1 - (1 - p[,2])^3
   p[,3] <- 1 - (1 - p[,3])^5
   for(i in 2:J) {
      pi[,i] <- pi[, i - 1]/p[, i - 1] * (1 - p[, i - 1]) * p[, i]
   }
   return(pi)
}
# .............................................................................


# pi function for capture-recapture design with 3 surveys
#    (introduced in AHM1 Section 7.8)
crPiFun <- function(p) {
   p1 <- p[,1]
   p2 <- p[,2]
   p3 <- p[,3]
   cbind("001" = (1 - p1) * (1 - p2) *      p3,
         "010" = (1 - p1) *      p2  * (1 - p3),
         "011" = (1 - p1) *      p2  *      p3,
         "100" =      p1  * (1 - p2) * (1 - p3),
         "101" =      p1  * (1 - p2) *      p3,
         "110" =      p1  *      p2  * (1 - p3),
         "111" =      p1  *      p2  *      p3)
}
# .............................................................................



# pi function for capture-recapture design with 3 surveys and behavioural response
#    (introduced in AHM1 Section 7.8.2)
crPiFun.Mb <- function(p) {
 pNaive <- p[,1]
 pWise <- p[,3]
 cbind("001" = (1 - pNaive) * (1 - pNaive) *      pNaive,
       "010" = (1 - pNaive) *      pNaive  * (1 - pWise),
       "011" = (1 - pNaive) *      pNaive  *      pWise,
       "100" =      pNaive  * (1 - pWise)  * (1 - pWise),
       "101" =      pNaive  * (1 - pWise)  *      pWise,
       "110" =      pNaive  *      pWise   * (1 - pWise),
       "111" =      pNaive  *      pWise   *      pWise)
}
# ..................................................................................


# Pi function for model with individual detection heterogeneity
#   (introduced in AHM1 Section 7.8.3)
MhPiFun <- function(p) {
  mu <- qlogis(p[,1])       # logit(p)
  sig <- exp(qlogis(p[1,2]))
  J <- ncol(p)
  M <- nrow(p)
  il <- matrix(NA, nrow=M, ncol=7)
  dimnames(il) <- list(NULL, c("001","010","011","100","101","110","111"))

for(i in 1:M) {
   il[i,1] <- integrate( function(x) {
     (1-plogis(mu[i]+x))*(1-plogis(mu[i]+x))*plogis(mu[i]+x)*dnorm(x,0,sig)
      }, lower=-Inf, upper=Inf, stop.on.error=FALSE)$value
   il[i,2] <- integrate( function(x) {
     (1-plogis(mu[i]+x))*plogis(mu[i]+x)*(1-plogis(mu[i]+x))*dnorm(x,0,sig)
      }, lower=-Inf, upper=Inf, stop.on.error=FALSE)$value
   il[i,3] <- integrate( function(x) {
      (1-plogis(mu[i]+x))*plogis(mu[i]+x)*plogis(mu[i]+x)*dnorm(x,0,sig)
      }, lower=-Inf, upper=Inf, stop.on.error=FALSE)$value
   il[i,4] <- integrate( function(x) {
      plogis(mu[i]+x)*(1-plogis(mu[i]+x))*(1-plogis(mu[i]+x))*dnorm(x,0,sig)
      }, lower=-Inf, upper=Inf, stop.on.error=FALSE)$value
   il[i,5] <- integrate( function(x) {
      plogis(mu[i]+x)*(1-plogis(mu[i]+x))*plogis(mu[i]+x)*dnorm(x,0,sig)
      }, lower=-Inf, upper=Inf, stop.on.error=FALSE)$value
   il[i,6] <- integrate( function(x) {
      plogis(mu[i]+x)*plogis(mu[i]+x)*(1-plogis(mu[i]+x))*dnorm(x,0,sig)
      }, lower=-Inf, upper=Inf, stop.on.error=FALSE)$value
   il[i,7] <- integrate( function(x) {
      plogis(mu[i]+x)*plogis(mu[i]+x)*plogis(mu[i]+x)*dnorm(x,0,sig)
      }, lower=-Inf, upper=Inf, stop.on.error=FALSE)$value
   }
return(il)
}
