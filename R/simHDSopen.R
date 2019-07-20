# Functions for the book Applied Hierarchical Modeling in Ecology (AHM)
# Marc Kery & Andy Royle, Academic Press, 2016.

# simHDSopen - AHM1 section 9.5.4.1 p499

# Function to generate open hierarchical distance sampling data
#   (introduced in AHM1 Section 9.5.4.1)

simHDSopen <-
function(type=c("line", "point"), nsites = 100, mean.lam = 2,
beta.lam = 0, mean.sig = 1, beta.sig = 0, B = 3, discard0=TRUE, nreps=2, phi=0.7, nyears=5, beta.trend = 0){
# Function simulates hierarchical distance sampling data under either a
#   line (type = "line") or a point (type = "point") transect protocol.
# Simulates a LIST OF LISTS of individual observations.
#   element [[i]][[j]] is the individual observations for replicate j of year i
# Function arguments:
#   nsites: Number of sites (spatial replication)
#   alpha.lam (= log(mean.lambda)), beta.lam: intercept and
#      slope of log-linear regression of expected lambda
#      on a habitat covariate
#   alpha.sig (= log(mean.sigma)), beta.sig: intercept and
#      slope of log-linear regression of scale parameter of
#      half-normal detection function on wind speed
#   B: strip half width
#
#   more things here
#
#  Note: for "point" the realized density is
#[(area of circle) /(area of  square)]*lambda

# Checks and fixes for input data -----------------------------
nsites <- round(nsites[1])
stopifNegative(mean.lam, allowZero=FALSE)
stopifNegative(B, allowZero=FALSE)
# --------------------------------------------

type <- match.arg(type)

parmvec <- c(mean.lam, beta.lam, mean.sig, beta.sig, phi, beta.trend)
names(parmvec)<- c("mean.lam", "beta.lam", "mean.sig", "beta.sig", "phi", "beta.trend")

# Make a covariate
habitat <- rnorm(nsites)                    # habitat covariate
# covariate "wind" is replicate level covariate
# Simulate abundance model (Poisson GLM for M)
M <-lambda <- matrix(NA, nrow=nsites, ncol=nyears)
Na <- wind <- array(NA,dim=c(nsites,nreps,nyears))
Na.real <- array(0, dim =c(nsites,nreps,nyears))
for(i in 1:nyears){
  lambda[,i] <- exp(log(mean.lam) + beta.lam*habitat + beta.trend*(i-nyears/2) )
  # density per "square"
  M[,i] <- rpois(nsites, lambda[,i])           # site-specific abundances
  Na[,,i] <- matrix(rbinom(nsites*nreps, M[,i],phi), nrow=nsites, byrow=FALSE)
  wind[,,i] <- runif(nsites*nreps, -2, 2)   # Wind covariate
}

# Detection probability model (site specific)
# this is now a 3-d array
sigma <- exp(log(mean.sig) + beta.sig*wind)

# Simulate observation model

## http://www.anderswallin.net/2009/05/uniform-random-points-in-a-circle-using-polar-coordinates/
outlist <-list()
for(yr in 1:nyears){

  list.yr <-list()
  for(rep in 1:nreps){
    data <- NULL
    for(i in 1:nsites){
      if(Na[i,rep,yr]==0){
        data <- rbind(data, c(i,NA,NA,NA,NA)) # save site, y=1, u1, u2, d
        next
      }
      if(type=="line"){
        # Simulation of distances, uniformly, for each ind. in pop.
        # note it piles up all M[i] guys on one side of the transect
        d <- runif(Na[i,rep,yr], 0, B)
        Na.real[i,rep,yr]<- sum(d<=B)
        p <- exp(-d *d / (2 * (sigma[i,rep,yr]^2)))
        # Determine if individuals are captured or not
        y <- rbinom(Na[i,rep,yr], 1, p)
        u1 <- u2 <- rep(NA, Na[i,rep,yr])   # coordinates (u,v)
        # Subset to "captured" individuals only
        d <- d[y==1]
        u1 <- u1[y==1]
        u2 <- u2[y==1]
        y <- y[y==1]
      }

      if(type=="point"){
        angle <- runif(Na[i,rep,yr], 0, 2*pi)
        dd<- B*sqrt(runif(Na[i,rep,yr],0,1))
        u1<- dd*cos(angle) + (B)
        u2<- dd*sin(angle) + (B)

        d <- sqrt((u1-B)^2 + (u2-B)^2)
        Na.real[i,rep,yr]<- sum(d<= B)
        p <- exp(-d *d / (2 * (sigma[i,rep,yr]^2)))
        # But we can only count individuals on a circle so we truncate p here
        pp <- ifelse(d < B, 1, 0) * p
        y <- rbinom(Na[i,rep,yr], 1, pp)  # Det./non-detection of each individual
        # Subset to "captured" individuals only
        u1 <- u1[y==1]
        u2 <- u2[y==1]
        d <- d[y==1]
        y <- y[y==1]
      }

      # Compile things into a matrix and insert NA if no individuals were
      # captured at site i. Coordinates (u,v) are not used here.
      if(sum(y) > 0) {
        data <- rbind(data, cbind(rep(i, sum(y)), y, u1, u2, d))
      } else {
        data <- rbind(data, c(i,NA,NA,NA,NA)) # make a row of missing data
      }
    } # end for(sites)
    colnames(data) <- c("site", "y", "u1", "u2", "d") # name 1st col "site"
    if(discard0)
      data <- data[!is.na(data[,2]),]
    list.yr[[rep]]<- data
  } # end for(rep)
  outlist[[yr]]<- list.yr
} # end for(year)
# Subset to sites at which individuals were captured. You may or may not
#  want to do this depending on how the model is formulated so be careful.
list(data=outlist, B=B, nsites=nsites, habitat=habitat, wind=wind, M.true= M, K=nreps,nyears=nyears,Na=Na, Na.real=Na.real,
  mean.lam=mean.lam, beta.lam=beta.lam, mean.sig=mean.sig, beta.sig=beta.sig, phi=phi, beta.trend=beta.trend, parms=parmvec )
}
