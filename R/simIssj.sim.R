# Functions for the book Applied Hierarchical Modeling in Ecology (AHM)
# Marc Kery & Andy Royle, Academic Press, 2016.

# issj.sim - AHM1 section 9.7.1 p517

# Function to simulate open distance sampling data for the Island Scrub Jays
#   (introduced in AHM1 Section 9.7.1)

issj.sim <-
function(B, db, lam, sigma, phi, gamma, npoints, nyrs, nbsize=-1.02){

stopifnot(min(db) == 0 && max(db) == B)
lam <- as.vector(lam)
sigma <- as.vector(sigma)
stopifnot(length(lam) == length(sigma))
nsites <- length(lam) # number of grid cells
stopifnot(phi >= 0 && phi <= 1)
stopifnot(gamma >= 0)
stopifnot(npoints <= nsites)

nD <- length(db) - 1 # Number of distance classes

Nsim <- matrix(0, nrow=nsites, ncol=nyrs) # simulated number of birds by site and year
#yr 1 as before
Nsim[,1]<-rnbinom (n=nsites, size=exp(nbsize), mu=lam) #generate individual counts per grid cell/point count circle
for(y in 2:nyrs){
  Nsim[,y]<-rbinom(nsites, Nsim[, y-1], phi) + rpois(nsites, Nsim[, y-1]*gamma)
}

#generate distance from hypothetical point count locations
# first, set prob for an individual to be in a 1m distance class from the center point
rc <- 1:B
ri <- (0:(B-1))
ar <- pi * (rc^2 - ri^2)
pcc <- ar/sum(ar)

NcList <- vector("list", nyrs)

for (y in 1:nyrs){
  NcList[[y]] <- matrix(nrow=0, ncol=2)
  for (j in 1:nsites){  # This is for all sites
    if (Nsim[j,y]==0) next
    junk <- rmultinom(1, Nsim[j,y], pcc) # count of birds in each 1m band
    tt <-rep( (which(junk!=0) - 0.5), (junk[which(junk!=0)]) ) # distance from centre
    Ndist <- cbind(rep(j,Nsim[j,y]), tt )
    NcList[[y]]<-rbind(NcList[[y]], Ndist)
  }
}

# for each sampling point generate detection data based on distance of individuals within a max of 300m
# and the detection model from the paper
cell<-sort(sample(1:nsites, npoints, replace=FALSE))

detList <- vector("list", nyrs)

for (y in 1:nyrs){
  for (j in cell) {
    dvec <- NcList[[y]][NcList[[y]][,1]==j, 2]
    if(length(dvec)==0) {
      det <- rep(0, nD)
    } else {
      pvec <- exp(-dvec^2/(2*(sigma[j]^2)))
      dets <- dvec[rbinom(length(dvec),1,pvec )==1]
      det <- table(cut(dets, db, include.lowest=TRUE))
    }
    detList[[y]]<-rbind(detList[[y]],det)
  }
}

# **npoints** x nyears matrix of total detections
y<-sapply(detList, rowSums)

# Pool all of the detection data into long vectors of distance category and site across all years
dclass<-site<-NULL
for (t in 1:nyrs){
  for (s in 1:npoints){
    if (y[s,t]==0) next

    ssi<-rep(cell[s], y[s,t])
    dc<-NULL
    for (k in 1:nD){
      if(detList[[t]][s,k]==0) next
      dd<-rep(k, detList[[t]][s,k])
      dc<-c(dc, dd)
    }
    dclass<-c(dclass, dc)
    site<-c(site, ssi)
  }
}


return(list(NcList=NcList, detList=detList, N=Nsim, cell=cell,
 y=y, dclass=dclass, site=site, nsites=nsites, lam=lam, phi=phi, gamma=gamma, sigma=sigma))

}
