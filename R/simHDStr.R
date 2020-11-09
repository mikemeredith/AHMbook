# Functions for the book Applied Hierarchical Modeling in Ecology (AHM)
# Marc Kery & Andy Royle, Academic Press, 2016.

# simHDStr - AHM1 section 9.3.2 p474


simHDStr <-
function(type = c("line", "point"), method=c("removal", "double"),
  nsites=200, lambda.group = 1, alpha0 = 0, alpha1 = 0,
  beta0 = 1, beta1 = 0.5, p.avail = 0.75, K = 3, p.double = c(0.4, 0.6), B = 3,
  discard0=FALSE, show.plot=TRUE){
# A general function for simulating hierarchical distance sampling (HDS) data
#   combined with a time-removal (with 3 removal periods) or
#   double-observer protocol, either for a line (type = "line") or
#   a point (type = "point") transect protocol
#   and with method = "removal" or method = "double".
#
# Other function arguments:
#     nsites: Number of sites (spatial replication)
#     lambda.group: Poisson mean of group size
#     alpha0, alpha1: intercept and slope of log-linear model relating sigma of
#        half-normal detection function to group size
#     beta0, beta1: intercept and slope of log-linear model relating the Poisson
#        mean of the number of groups per unit area to habitat
#     p.avail: overall availability probability (phi in text)
#int.avail: time interval-specific availability probability
#     K: number of removal periods (of equal length)
#     p.double: detection probability for first and second observer
#     B: strip half width
#     discard0: whether to discard or keep the data from sites with nobody detected

# Checks and fixes for input data -----------------------------
nsites <- round(nsites[1])
stopifNegative(lambda.group, allowZero=FALSE)
stopifnotProbability(p.avail)
K <- round(K)[1]
stopifnotProbability(p.double)
stopifNegative(B, allowZero=FALSE)
# --------------------------------------------

type <- match.arg(type)
method <- match.arg(method)

# Create covariate
habitat <- rnorm(nsites)              # Simulated continuous covariate

# Simulate superpopulation abundance model for groups (Poisson GLM for M)
lambda <- exp(beta0 + beta1*habitat)  # Density of groups per "square"
M <- rpois(nsites, lambda)            # site-specific number of groups
M.true <- M                           # for point: inside of B

data <- NULL
for(i in 1:nsites){
  if(M[i]==0){
    data <- rbind(data,c(i,NA,NA,NA,NA,NA,NA)) # save site, y=1, u1, u2, d, gs, tint
  next
  }

  # Simulation for line transect data
  if(type=="line"){
    # Simulation of distances, uniformly, for each individual in the population
    # note it piles up all N[i] guys on one side of the transect
    d <- runif(M[i], 0, B)
    gs <- rpois(M[i], lambda.group) + 1    # Observable group size
    sigma.vec<- exp(alpha0 + alpha1*(gs-1) ) # subtract 1 for interpretation
    # Detection probability for each group
    p <- exp(-d*d/(2*(sigma.vec^2)))
    # Time-removal protocol
    if(method=="removal"){
      int.avail <- 1 - (1-p.avail)^(1/K)
      rem.probs <- c(int.avail, ((1-int.avail)^(1:(K-1)))*int.avail)
      mn.probs <- c(rem.probs, 1-sum(rem.probs))
      # was this (think wrong):   aux <- sample(1:(K+1), N[i], replace=TRUE, prob=mn.probs)
      aux <- sample(1:(K+1), M[i], replace=TRUE, prob=mn.probs)
      aux[aux==(K+1)] <- 0
    }
    # Double-observer protocol
    if(method=="double"){
      rem.probs <- c(p.double[1]*(1-p.double[2]), (1-p.double[1]) * p.double[2],
        p.double[1]*p.double[2])
      mn.probs <- c(rem.probs, 1-sum(rem.probs))
      aux <- sample(1:4, M[i], replace=TRUE, prob=mn.probs)
      aux[aux==4]<- 0
    }
    newp <-  p * as.numeric(aux!=0)
    navail <- sum(aux!=0)

    if(navail==0){
      data <- rbind(data,c(i,NA,NA,NA,NA,NA,NA)) # save site, y=1, u1, u2, d
      next
    }

    # generate count of birds based on combined probability of detection
    y <- rbinom(M[i], 1, newp)
    # Subset to "captured" individuals only
    u1 <- u2 <- rep(NA, sum(y))  ;  d <- d[y==1]  ;  gs <- gs[y==1]
    aux <- aux[y==1]  ;  y <- y[ y==1]
  }

  # Simulation for point transect data
  if(type=="point"){
     # Simulation of data on a circle of radius B
     angle <- runif(M[i], 0, 2*pi)
     r2 <- runif(M[i], 0, 1)
     r<-  B*sqrt(r2)
     u1<-  r*cos(angle) + B
     u2<-  r*sin(angle) + B

    d <- sqrt((u1 - B)^2 + (u2-B)^2)  ## d == r ! Cruft...
    M.true[i] <- sum(d<= B)    # Population size inside of count circle
    gs <- rpois(M[i], lambda.group) + 1
    sigma.vec <- exp(alpha0 + alpha1*(gs-1))
    # But we can only count individuals on a circle so we truncate p here
    p <- ifelse(d<(B),1,0)*exp(-d*d/(2*(sigma.vec^2)))

    # Time-removal protocol
    if(method=="removal"){
      int.avail <- 1 - (1-p.avail)^(1/K)
      rem.probs <- c(int.avail, ((1-int.avail)^(1:(K-1)))*int.avail)
      mn.probs <- c(rem.probs, 1-sum(rem.probs))
      aux <- sample(1:(K+1), M[i], replace=TRUE, prob=mn.probs)
      aux[aux==(K+1)] <- 0
    }
    # Double-observer protocol
    if(method=="double"){
      rem.probs <- c(p.double[1]*(1-p.double[2]), (1-p.double[1])*p.double[2],
        p.double[1]*p.double[2])
      mn.probs <- c(rem.probs, 1-sum(rem.probs))
      aux <- sample(1:(K+1), M[i], replace=TRUE, prob=mn.probs)
      aux[aux==(K+1)]<- 0
    }

    newp <-  p * as.numeric(aux!=0)
    navail <- sum(aux!=0)

    if(navail==0){
      data <- rbind(data,c(i,NA,NA,NA,NA,NA,NA)) # save site, y=1, u1, u2, d
      next
    }

    # generate count of birds based on combined probability of detection
    y <- rbinom(M[i], 1, newp)
    # Subset to "captured" individuals only
    u1 <- u1[y==1]  ;  u2 <- u2[y==1]  ;  d <- d[y==1]  ;  gs <- gs[ y==1]
    aux <- aux[y==1]  ;  y <- y[ y==1]
  }

  # Now compile things into a matrix and insert NA if no individuals were
  # captured at site i. Coordinates (u,v) are not used here.
  if(sum(y)>0){
    data <- rbind(data, cbind(rep(i, sum(y)), y, u1, u2, d, gs, aux))
  } else {
    data <- rbind(data, c(i,NA,NA,NA,NA,NA,NA)) # make a row of missing data
  }
} # end of for loop
colnames(data)[1] <- "site"
# Subset to sites at which individuals were captured. You may or may not
#  do this depending on how the model is formulated so be careful.
if(discard0)
  data <- data[!is.na(data[,2]),]
parmvec <- c(alpha0, alpha1, beta0, beta1, p.avail, p.double)
names(parmvec) <- c("alpha0", "alpha1", "beta0", "beta1", "p.avail",
  "p.double1", "p.double2")

# Visualisation
if(show.plot) {
  if(type=="line"){       # For line transect
    op <- par(mfrow = c(1, 3)) ; on.exit(par(op))
    tryPlot <- try( {
      hist(data[,"d"], col = "lightblue", breaks = 20,
          main = "Frequency of distances to groups", xlab = "Distance")
      ttt <- table(data[,1])
      n <- rep(0, nsites)
      n[as.numeric(rownames(ttt))] <- ttt
      plot(habitat, n, main = "Observed group counts (n) vs. habitat", frame = FALSE)
      plot(table(data[,"gs"]), main = "Observed group sizes", ylab = "Frequency",
          frame = FALSE)
    }, silent = TRUE)
    if(inherits(tryPlot, "try-error"))
      tryPlotError(tryPlot)
  }

  if(type=="point"){       # For point transect
    op <- par(mfrow = c(2,2)) ; on.exit(par(op))
    tryPlot <- try( {
      plot(data[,"u1"], data[,"u2"], pch = 16,
          main = "Located groups in point transects", xlim = c(0, 2*B),
          ylim = c(0, 2*B), col = data[,1], asp = 1)
      points(B, B, pch = "+", cex = 3)
      # library(plotrix)
      draw.circle(B, B, B)
      hist(data[,"d"], col = "lightblue", breaks = 20, main =
        "Frequency of distances to groups", xlab = "Distance")
      ttt <- table(data[,1])
      n <- rep(0, nsites)
      n[as.numeric(rownames(ttt))] <- ttt
      plot(habitat, n, main = "Observed group counts (n) vs. habitat", frame = FALSE)
      plot(table(data[,"gs"]), main = "Observed group sizes", ylab = "Frequency",
          frame = FALSE)
    }, silent = TRUE)
    if(inherits(tryPlot, "try-error"))
      tryPlotError(tryPlot)
  }
}

# Output
list(type = type, method = method, nsites = nsites, lambda.group = lambda.group, alpha0 = alpha0, alpha1 = alpha1, beta0 = beta0, beta1 = beta1, p.avail = p.avail, p.double = p.double, K = K, B = B, data=data, habitat=habitat, M = M, M.true = M.true, parms = parmvec)
}
