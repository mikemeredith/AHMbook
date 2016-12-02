\name{issj.sim}
\alias{issj.sim}
\encoding{UTF-8}
\title{
Simulate open distance sampling data for the Island Scrub Jays
}
\description{
Function to simulate open distance sampling data for the Island Scrub Jays
}
\usage{
issj.sim(B = 300, db = db, lam = lamnew, sigma = sigma, phi = phi, gamma = gamma, npoints = npoints, nyrs = nyrs, nbsize = -1.02)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{B}{
%%     ~~Describe \code{B} here~~
}
  \item{db}{
%%     ~~Describe \code{db} here~~
}
  \item{lam}{
%%     ~~Describe \code{lam} here~~
}
  \item{sigma}{
%%     ~~Describe \code{sigma} here~~
}
  \item{phi}{
%%     ~~Describe \code{phi} here~~
}
  \item{gamma}{
%%     ~~Describe \code{gamma} here~~
}
  \item{npoints}{
%%     ~~Describe \code{npoints} here~~
}
  \item{nyrs}{
%%     ~~Describe \code{nyrs} here~~
}
  \item{nbsize}{
%%     ~~Describe \code{nbsize} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (B = 300, db = db, lam = lamnew, sigma = sigma, phi = phi, 
    gamma = gamma, npoints = npoints, nyrs = nyrs, nbsize = -1.02) 
{
    parms <- c(lam = lam, phi = phi, gamma = gamma, sigma = sigma)
    J = dim(lam)[1]
    cell <- sort(sample(1:J, npoints, replace = FALSE))
    Nsim <- matrix(0, nrow = npoints, ncol = nyrs)
    Nsim[, 1] <- rnbinom(n = npoints, size = exp(nbsize), mu = lam[cell, 
        1])
    for (y in 2:nyrs) {
        Nsim[, y] <- rbinom(npoints, Nsim[, y - 1], phi) + rpois(npoints, 
            Nsim[, y - 1] * gamma)
    }
    rc <- 1:B
    ri <- (0:(B - 1))
    ar <- NULL
    for (i in 1:B) {
        ar[i] <- pi * (rc[i]^2 - ri[i]^2)
    }
    pcc <- ar/sum(ar)
    Nc <- matrix(nrow = 0, ncol = 2)
    NcList <- list(Nc, Nc, Nc, Nc, Nc, Nc)
    for (y in 1:nyrs) {
        for (j in 1:J) {
            if (Nsim[j, y] == 0) 
                next
            junk <- rmultinom(1, Nsim[j, y], pcc)
            tt <- rep((which(junk != 0) - 0.5), (junk[which(junk != 
                0)]))
            Ndist <- cbind(rep(j, Nsim[j, y]), tt)
            NcList[[y]] <- rbind(NcList[[y]], Ndist)
        }
    }
    detec <- NULL
    detList <- list(detec, detec, detec, detec, detec, detec)
    for (y in 1:nyrs) {
        for (j in 1:J) {
            dvec <- NcList[[y]][which(NcList[[y]][, 1] == j), 
                2]
            if (length(dvec) == 0) {
                det <- c(rep(0, length(db) - 1))
            }
            else {
                pvec <- exp(-dvec * dvec/(2 * (sigma[j]^2)))
                dets <- dvec[rbinom(length(dvec), 1, pvec) == 
                  1]
                det <- table(cut(dets, db, include.lowest = T))
            }
            detList[[y]] <- rbind(detList[[y]], det)
        }
    }
    y <- sapply(detList, rowSums)
    dclass <- site <- NULL
    for (t in 1:nyrs) {
        for (s in 1:npoints) {
            if (y[s, t] == 0) 
                next
            ssi <- rep(s, y[s, t])
            dc <- NULL
            for (k in 1:nD) {
                if (detList[[t]][s, k] == 0) 
                  next
                dd <- rep(k, detList[[t]][s, k])
                dc <- c(dc, dd)
            }
            dclass <- c(dclass, dc)
            site <- c(site, ssi)
        }
    }
    return(list(NcList = NcList, detList = detList, N = Nsim, 
        cell = cell, parms = parms, y = y, dclass = dclass, site = site))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
