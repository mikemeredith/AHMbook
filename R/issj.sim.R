# Functions for the book Applied Hierarchical Modeling in Ecology (AHM)
# Marc Kery & Andy Royle, Academic Press, 2016.

# issj.sim - section 9.7.1 p517

# Function to simulate open distance sampling data for the Island Scrub Jays
#   (introduced in Section 9.7.1)

issj.sim <-
function(B, db, lam, sigma, phi, gamma, npoints, nyrs, nbsize=-1.02){

nD <- length(db) - 1 # Number of distance classes

parms<- c(lam=lam,phi=phi, gamma=gamma, sigma=sigma)

J=dim(lam)[1] #number of grid cells

cell<-sort(sample(1:J, npoints, replace=FALSE))

Nsim<-matrix(0,nrow=npoints, ncol=nyrs)

#yr 1 as before
Nsim[,1]<-rnbinom (n=npoints, size=exp(nbsize), mu=lam[cell,1]) #generate individual countes per grid cell/point count circle
for(y in 2:nyrs){
  Nsim[,y]<-rbinom(npoints,Nsim[,y-1], phi) + rpois(npoints,Nsim[,y-1]*gamma)
}
#generate distance from hypothetical point count locations
#first set prob for an individual to be in a 1-m distance class from the center point

rc<-1:B
ri<-(0:(B-1))
ar<-NULL
for (i in 1:B) {
   ar[i]<-pi * (rc[i]^2 - ri[i]^2)
}
pcc<-ar/sum(ar)

Nc<-matrix(nrow=0, ncol=2)
NcList<-list(Nc,Nc,Nc,Nc,Nc,Nc)

for (y in 1:nyrs){
for (j in 1:J){
  if (Nsim[j,y]==0) next
  junk <- rmultinom(1, Nsim[j,y], pcc)
  tt <-rep( (which(junk!=0) - 0.5), (junk[which(junk!=0)]) )
  Ndist<-cbind(rep(j,Nsim[j,y]),tt )
  NcList[[y]]<-rbind(NcList[[y]], Ndist)
}}

# for each sampling point generate detection data based on distance of individuals within a max of 300m
# and the detection model from the paper

detec<-NULL
detList<-list(detec, detec,detec,detec,detec,detec)

for (y in 1:nyrs){
for (j in 1:J) {
  dvec<-NcList[[y]][ which(NcList[[y]][,1]==j)  ,2 ]
  if(length(dvec)==0) {
   det<-c(rep(0,length(db)-1))
   } else {
   pvec<-exp(-dvec*dvec/(2*(sigma[j]^2)))
   dets<-dvec[rbinom(length(dvec),1,pvec )==1]
   det<-table(cut(dets, db, include.lowest=T))
  }
  detList[[y]]<-rbind(detList[[y]],det)
 }
}

# nsites x nyears matrix of total detections
y<-sapply(detList, rowSums)

# Pool all of the detection data into long vectors of distance category and site across all years
dclass<-site<-NULL
for (t in 1:nyrs){
 for (s in 1:npoints){
    if (y[s,t]==0) next

    ssi<-rep(s, y[s,t])
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


return(list(NcList=NcList, detList=detList, N=Nsim, cell=cell, parms=parms, y=y, dclass=dclass, site=site))

}
