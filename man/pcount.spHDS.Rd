\name{pcount.spHDS}
\alias{pcount.spHDS}
\encoding{UTF-8}

\title{
Fit spatial hierarchical distance sampling model
}
\description{
Function fits spatial hierarchical distance sampling model unmarked model fitting function introduced in Section 9.8.4.
}
\usage{
pcount.spHDS(formula, data, K, mixture = c("P", "NB", "ZIP"), starts, method = "BFGS", se = TRUE, engine = c("C", "R"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
%%     ~~Describe \code{formula} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{K}{
%%     ~~Describe \code{K} here~~
}
  \item{mixture}{
%%     ~~Describe \code{mixture} here~~
}
  \item{starts}{
%%     ~~Describe \code{starts} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{se}{
%%     ~~Describe \code{se} here~~
}
  \item{engine}{
%%     ~~Describe \code{engine} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
function (formula, data, K, mixture = c("P", "NB", "ZIP"), starts, 
    method = "BFGS", se = TRUE, engine = c("C", "R"), ...) 
{
    mixture <- match.arg(mixture, c("P", "NB", "ZIP"))
    if (!is(data, "unmarkedFramePCount")) 
        stop("Data is not an unmarkedFramePCount object.")
    engine <- match.arg(engine, c("C", "R"))
    if (identical(mixture, "ZIP") & identical(engine, "R")) 
        stop("ZIP mixture not available for engine='R'")
    designMats <- unmarked:::getDesign(data, formula)
    X <- designMats$X
    V <- designMats$V
    y <- designMats$y
    X.offset <- designMats$X.offset
    V.offset <- designMats$V.offset
    if (is.null(X.offset)) 
        X.offset <- rep(0, nrow(X))
    if (is.null(V.offset)) 
        V.offset <- rep(0, nrow(V))
    NAmat <- is.na(y)
    J <- ncol(y)
    M <- nrow(y)
    lamParms <- colnames(X)
    detParms <- colnames(V)
    nDP <- ncol(V)
    nAP <- ncol(X)
    if (missing(K)) {
        K <- max(y, na.rm = TRUE) + 100
        warning("K was not specified and was set to ", K, ".")
    }
    if (K <= max(y, na.rm = TRUE)) 
        stop("specified K is too small. Try a value larger than any observation")
    k <- 0:K
    lk <- K + 1
    M <- nrow(y)
    J <- ncol(y)
    k.ik <- rep(k, M)
    k.ijk <- rep(k, M * J)
    nP <- nAP + nDP + (mixture != "P")
    if (!missing(starts) && length(starts) != nP) 
        stop(paste("The number of starting values should be", 
            nP))
    if (identical(engine, "R")) {
        y.ij <- as.numeric(t(y))
        y.ijk <- rep(y.ij, each = K + 1)
        navec <- is.na(y.ijk)
        ijk <- expand.grid(k = 0:K, j = 1:J, i = 1:M)
        ijk.to.ikj <- with(ijk, order(i, k, j))
        nll <- function(parms) {
            theta.i <- exp(X \%*\% parms[1:nAP] + X.offset)
            p.ij <- 2 * plogis(V \%*\% parms[(nAP + 1):(nAP + nDP)] + 
                V.offset)
            theta.ik <- rep(theta.i, each = K + 1)
            p.ijk <- rep(p.ij, each = K + 1)
            bin.ijk <- dbinom(y.ijk, k.ijk, p.ijk)
            bin.ijk[which(is.na(bin.ijk))] <- 1
            bin.ik.mat <- matrix(bin.ijk[ijk.to.ikj], M * (K + 
                1), J, byrow = TRUE)
            g.ik <- rowProds(bin.ik.mat)
            if (identical(mixture, "P")) {
                f.ik <- dpois(k.ik, theta.ik)
            }
            else if (identical(mixture, "NB")) {
                f.ik <- dnbinom(k.ik, mu = theta.ik, size = exp(parms[nP]))
            }
            dens.i.mat <- matrix(f.ik * g.ik, M, K + 1, byrow = TRUE)
            dens.i <- rowSums(dens.i.mat)
            -sum(log(dens.i))
        }
    }
    else {
        nll <- function(parms) {
            beta.lam <- parms[1:nAP]
            beta.p <- parms[(nAP + 1):(nAP + nDP)]
            log.alpha <- 1
            if (mixture \%in\% c("NB", "ZIP")) 
                log.alpha <- parms[nP]
            .Call("nll_pcount", y, X, V, beta.lam, beta.p, log.alpha, 
                X.offset, V.offset, NAmat, lk, mixture, PACKAGE = "unmarked")
        }
    }
    if (missing(starts)) 
        starts <- rep(0, nP)
    fm <- optim(starts, nll, method = method, hessian = se, ...)
    opt <- fm
    ests <- fm$par
    nbParm <- switch(mixture, NB = "alpha", ZIP = "psi", P = character(0))
    names(ests) <- c(lamParms, detParms, nbParm)
    if (se) {
        tryCatch(covMat <- solve(fm$hessian), error = function(x) stop(simpleError("Hessian is singular.  Try using fewer covariates.")))
    }
    else {
        covMat <- matrix(NA, nP, nP)
    }
    fmAIC <- 2 * fm$value + 2 * nP
    stateName <- "Abundance"
    stateEstimates <- unmarked:::unmarkedEstimate(name = stateName, 
        short.name = "lam", estimates = ests[1:nAP], covMat = as.matrix(covMat[1:nAP, 
            1:nAP]), invlink = "exp", invlinkGrad = "exp")
    detEstimates <- unmarked:::unmarkedEstimate(name = "Detection", 
        short.name = "p", estimates = ests[(nAP + 1):(nAP + nDP)], 
        covMat = as.matrix(covMat[(nAP + 1):(nAP + nDP), (nAP + 
            1):(nAP + nDP)]), invlink = "logistic", invlinkGrad = "logistic.grad")
    estimateList <- unmarked:::unmarkedEstimateList(list(state = stateEstimates, 
        det = detEstimates))
    if (identical(mixture, "NB")) {
        estimateList@estimates$alpha <- unmarked:::unmarkedEstimate(name = "Dispersion", 
            short.name = "alpha", estimates = ests[nP], covMat = as.matrix(covMat[nP, 
                nP]), invlink = "exp", invlinkGrad = "exp")
    }
    if (identical(mixture, "ZIP")) {
        estimateList@estimates$psi <- unmarked:::unmarkedEstimate(name = "Zero-inflation", 
            short.name = "psi", estimates = ests[nP], covMat = as.matrix(covMat[nP, 
                nP]), invlink = "logistic", invlinkGrad = "logistic.grad")
    }
    umfit <- new("unmarkedFitPCount", fitType = "pcount", call = match.call(), 
        formula = formula, data = data, sitesRemoved = designMats$removed.sites, 
        estimates = estimateList, AIC = fmAIC, opt = opt, negLogLike = fm$value, 
        nllFun = nll, K = K, mixture = mixture)
    return(umfit)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
