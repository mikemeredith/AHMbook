\name{simHDStr}
\alias{simHDStr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
simulate HDS time-removal data%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
simHDStr(type = "line", method = "removal", nsites = 200, lambda.group = 1, alpha0 = 0, alpha1 = 0, beta0 = 1, beta1 = 0.5, p.avail = 0.75, K = 3, p.double = c(0.4, 0.6), B = 3, discard0 = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{type}{
%%     ~~Describe \code{type} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{nsites}{
%%     ~~Describe \code{nsites} here~~
}
  \item{lambda.group}{
%%     ~~Describe \code{lambda.group} here~~
}
  \item{alpha0}{
%%     ~~Describe \code{alpha0} here~~
}
  \item{alpha1}{
%%     ~~Describe \code{alpha1} here~~
}
  \item{beta0}{
%%     ~~Describe \code{beta0} here~~
}
  \item{beta1}{
%%     ~~Describe \code{beta1} here~~
}
  \item{p.avail}{
%%     ~~Describe \code{p.avail} here~~
}
  \item{K}{
%%     ~~Describe \code{K} here~~
}
  \item{p.double}{
%%     ~~Describe \code{p.double} here~~
}
  \item{B}{
%%     ~~Describe \code{B} here~~
}
  \item{discard0}{
%%     ~~Describe \code{discard0} here~~
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
function (type = "line", method = "removal", nsites = 200, lambda.group = 1, 
    alpha0 = 0, alpha1 = 0, beta0 = 1, beta1 = 0.5, p.avail = 0.75, 
    K = 3, p.double = c(0.4, 0.6), B = 3, discard0 = FALSE) 
{
    habitat <- rnorm(nsites)
    lambda <- exp(beta0 + beta1 * habitat)
    M <- rpois(nsites, lambda)
    M.true <- M
    data <- NULL
    for (i in 1:nsites) {
        if (M[i] == 0) {
            data <- rbind(data, c(i, NA, NA, NA, NA, NA, NA))
            next
        }
        if (type == "line") {
            d <- runif(M[i], 0, B)
            gs <- rpois(M[i], lambda.group) + 1
            sigma.vec <- exp(alpha0 + alpha1 * (gs - 1))
            p <- exp(-d * d/(2 * (sigma.vec^2)))
            if (method == "removal") {
                int.avail <- 1 - (1 - p.avail)^(1/K)
                rem.probs <- c(int.avail, ((1 - int.avail)^(1:(K - 
                  1))) * int.avail)
                mn.probs <- c(rem.probs, 1 - sum(rem.probs))
                aux <- sample(1:(K + 1), M[i], replace = TRUE, 
                  prob = mn.probs)
                aux[aux == (K + 1)] <- 0
            }
            if (method == "double") {
                rem.probs <- c(p.double[1] * (1 - p.double[2]), 
                  (1 - p.double[1]) * p.double[2], p.double[1] * 
                    p.double[2])
                mn.probs <- c(rem.probs, 1 - sum(rem.probs))
                aux <- sample(1:4, M[i], replace = TRUE, prob = mn.probs)
                aux[aux == 4] <- 0
            }
            newp <- p * as.numeric(aux != 0)
            navail <- sum(aux != 0)
            if (navail == 0) {
                data <- rbind(data, c(i, NA, NA, NA, NA, NA, 
                  NA))
                next
            }
            y <- rbinom(M[i], 1, newp)
            u1 <- u2 <- rep(NA, sum(y))
            d <- d[y == 1]
            gs <- gs[y == 1]
            aux <- aux[y == 1]
            y <- y[y == 1]
        }
        if (type == "point") {
            angle <- runif(M[i], 0, 360)
            r2 <- runif(M[i], 0, 1)
            r <- B * sqrt(r2)
            u1 <- r * cos(angle) + B
            u2 <- r * sin(angle) + B
            d <- sqrt((u1 - B)^2 + (u2 - B)^2)
            M.true[i] <- sum(d <= B)
            gs <- rpois(M[i], lambda.group) + 1
            sigma.vec <- exp(alpha0 + alpha1 * (gs - 1))
            p <- ifelse(d < (B), 1, 0) * exp(-d * d/(2 * (sigma.vec^2)))
            if (method == "removal") {
                int.avail <- 1 - (1 - p.avail)^(1/K)
                rem.probs <- c(int.avail, ((1 - int.avail)^(1:(K - 
                  1))) * int.avail)
                mn.probs <- c(rem.probs, 1 - sum(rem.probs))
                aux <- sample(1:(K + 1), M[i], replace = TRUE, 
                  prob = mn.probs)
                aux[aux == (K + 1)] <- 0
            }
            if (method == "double") {
                rem.probs <- c(p.double[1] * (1 - p.double[2]), 
                  (1 - p.double[1]) * p.double[2], p.double[1] * 
                    p.double[2])
                mn.probs <- c(rem.probs, 1 - sum(rem.probs))
                aux <- sample(1:(K + 1), M[i], replace = TRUE, 
                  prob = mn.probs)
                aux[aux == (K + 1)] <- 0
            }
            newp <- p * as.numeric(aux != 0)
            navail <- sum(aux != 0)
            if (navail == 0) {
                data <- rbind(data, c(i, NA, NA, NA, NA, NA, 
                  NA))
                next
            }
            y <- rbinom(M[i], 1, newp)
            u1 <- u1[y == 1]
            u2 <- u2[y == 1]
            d <- d[y == 1]
            gs <- gs[y == 1]
            aux <- aux[y == 1]
            y <- y[y == 1]
        }
        if (sum(y) > 0) {
            data <- rbind(data, cbind(rep(i, sum(y)), y, u1, 
                u2, d, gs, aux))
        }
        else data <- rbind(data, c(i, NA, NA, NA, NA, NA, NA))
    }
    if (discard0) 
        data <- data[!is.na(data[, 2]), ]
    parmvec <- c(alpha0, alpha1, beta0, beta1, p.avail, p.double)
    names(parmvec) <- c("alpha0", "alpha1", "beta0", "beta1", 
        "p.avail", "p.double1", "p.double2")
    if (type == "line") {
        par(mfrow = c(1, 3))
        hist(data[, "d"], col = "lightblue", breaks = 20, main = "Frequency of distances to groups", 
            xlab = "Distance")
        ttt <- table(data[, 1])
        n <- rep(0, nsites)
        n[as.numeric(rownames(ttt))] <- ttt
        plot(habitat, n, main = "Observed group counts (n) vs. habitat", 
            frame = F)
        plot(table(data[, "gs"]), main = "Observed group sizes", 
            ylab = "Frequency", frame = F)
    }
    if (type == "point") {
        par(mfrow = c(2, 2))
        plot(data[, "u1"], data[, "u2"], pch = 16, main = "Located groups in point transects", 
            xlim = c(0, 2 * B), ylim = c(0, 2 * B), col = data[, 
                1], asp = 1)
        points(B, B, pch = "+", cex = 3)
        library(plotrix)
        draw.circle(B, B, B)
        hist(data[, "d"], col = "lightblue", breaks = 20, main = "Frequency of distances to groups", 
            xlab = "Distance")
        ttt <- table(data[, 1])
        n <- rep(0, nsites)
        n[as.numeric(rownames(ttt))] <- ttt
        plot(habitat, n, main = "Observed group counts (n) vs. habitat", 
            frame = F)
        plot(table(data[, "gs"]), main = "Observed group sizes", 
            ylab = "Frequency", frame = F)
    }
    list(type = type, method = method, nsites = nsites, lambda.group = lambda.group, 
        alpha0 = alpha0, alpha1 = alpha1, beta0 = beta0, beta1 = beta1, 
        p.avail = p.avail, p.double = p.double, K = K, B = B, 
        data = data, habitat = habitat, M = M, M.true = M.true, 
        parms = parmvec)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
