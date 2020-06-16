
# Code from Andy, 29 Dec 2016

# Mike reorganised the plotting commands, adding devAskNewPage and show.plots.
# Mike moved beta1 and npix to the arguments.

sim.spatialHDS <-
function(lam0 = 4 , sigma= 1.5,  B=3, nsites=100, beta1 = 1, npix = 20, show.plots = 3){

# Function simulates coordinates of individuals on a square
# Square is [0,2B] x[0,2B], with a count location on the point (B,B)
# lam0: expected population size in the square
# sigma: scale of half-normal detection function
# B: circle radius

# Checks and fixes for input data -----------------------------
stopifNegative(lam0, allowZero=FALSE)
stopifNegative(sigma, allowZero=FALSE)
stopifNegative(B, allowZero=FALSE)
nsites <- round(nsites[1])
npix <- round(npix[1])
# --------------------------------------------

if(show.plots > 0) {
  oldpar <- par(mar=c(3,3,3,6), "mfrow")
  oldAsk <- devAskNewPage(ask = dev.interactive(orNone = TRUE))
  on.exit({par(oldpar) ; devAskNewPage(oldAsk)})
}

# npix<- 20
data<- NULL
beta0<- log(lam0/(npix*npix))
# beta1<- 1


Z<- matrix(NA,nrow=npix*npix, ncol=nsites)

delta<- (2*B-0)/npix
grx<- seq(delta/2, 2*B - delta/2, delta)
gr<- expand.grid(grx,grx)
V<- exp(-e2dist(gr,gr)/1)
N<- rep(NA,nsites)

for(s in 1:nsites){
  z<- t(chol(V))%*%rnorm( npix^2 )
  Z[,s]<- z

  # Note Poisson assumption which means in each pixel is also Poisson
  N[s]<- rpois(1, sum(exp( beta0 + beta1*Z[,s])))

  probs<- exp(beta1*Z[,s])/sum(exp(beta1*Z[,s]))
  pixel.id<- sample(1:(npix^2), N[s], replace=TRUE, prob=probs)
  # could simulate ranomdly within the pixel but it won't matter
  u1<- gr[pixel.id,1]
  u2<- gr[pixel.id,2]

  d <- sqrt((u1 - B)^2 + (u2-B)^2) # distance to center point of square

  p<- exp(-d*d/(2*sigma*sigma))

  # Now we decide whether each individual is detected or not
  y <- rbinom(N[s], 1, p)

  if(s <= show.plots) {
    tryPlot <- try( {
      img <- rasterFromXYZ(cbind(gr,z))
      image(img, col=topo.colors(10))
      #draw.circle(3,3,B)
      image_scale(z,col=topo.colors(10))
      points(u1,u2,pch=16,col='black')

      # points(u1[d<= B], u2[d<= B], pch = 16, col = "black")
      points(u1[y==1], u2[y==1], pch = 16, col = "red")
      points(B, B,   ,pch = "+", cex = 3)
      # draw.circle(3, 3,   B)
    }, silent = TRUE)
    if(inherits(tryPlot, "try-error")) {
      show.plots <- 0 # stop further plotting attempts
      tryPlotError(tryPlot)
    }
  }

  if(sum(y)>0) {
    data<- rbind(data, cbind(rep(s,length(u1)),u1=u1,u2=u2,d=d,y=y))
  } else {
    data<- rbind(data, c(s, NA, NA, NA, NA))
  }
}

dimnames(data)<-list(NULL,c("site","u1","u2","d","y"))

return(list(data=data, B=B, Habitat=Z, grid=gr,N=N,nsites=nsites))
}
