\name{simComm}
\alias{simComm}
\encoding{UTF-8}
\title{
Simulate community occupancy or community abundance data
}
\description{
Function to simulate community occupancy or community abundance data with random species effects for psi/lambda and p (both including effects of one covariate, 'habitat' for psi/lambda and 'wind speed' for p) (introduced in Section 11.2)
}
\usage{
simComm(type = "det/nondet", nsite = 30, nrep = 3, nspec = 100, 
  mean.psi = 0.25, sig.lpsi = 1, mu.beta.lpsi = 0, sig.beta.lpsi = 0,
  mean.lambda = 2, sig.loglam = 1, mu.beta.loglam = 1, sig.beta.loglam = 1,
  mean.p = 0.25, sig.lp = 1, mu.beta.lp = 0, sig.beta.lp = 0, show.plot = TRUE)
}
\arguments{
  \item{type}{
"det/nondet" or "counts"; choose whether you want to simulate detection/nondetection data or count data
}
  \item{nsite}{
number of sites
}
  \item{nrep}{
number of replicate samples (occasions or repeated measurements)
}
  \item{nspec}{
total number of species in the area that is sampled by these sites (regional species pool)
}
  \item{mean.psi}{
community mean of occupancy probability over all species in community (probability scale)
}
  \item{sig.lpsi}{
 community standard deviation of logit(occupancy probability intercept)
}
  \item{mu.beta.lpsi}{
 community mean of the effects of 'habitat' covariate on occupancy probability
}
  \item{sig.beta.lpsi}{
community standard deviation of the effects of 'habitat' covariate on occupancy probability
}
  \item{mean.lambda}{
community mean of expected abundance over all species in superpopulation
}
  \item{sig.loglam}{
community standard deviation of log(lambda intercept)
}
  \item{mu.beta.loglam}{
 community mean of the effects of 'habitat' covariate on log(lambda)
}
  \item{sig.beta.loglam}{
community standard deviation of the effects of 'habitat' covariate on expected abundance
}
  \item{mean.p}{
community mean of detection probability over all species in superpopulation (probability scale)
}
  \item{sig.lp}{
community standard deviation of logit(detection probability intercept)
}
  \item{mu.beta.lp}{
community mean of the effects of 'wind' covariate on detection probability
}
  \item{sig.beta.lp}{
community standard deviation of the effects of 'wind'covariate on detection probability
}
  \item{show.plot}{
choose whether to show plots or not. Set to FALSE when using function in simulations.
}
}
\details{
Function simulates data from repeated sampling of a metacommunity (or spatially structured community) according the model of Dorazio & Royle (JASA, 2005) for type = "det/nondet" (this is the default) or under the model of Yamaura et al. (2012) for type = "counts".

Occupancy probability (psi) or expected abundance (lambda) can be made dependent on a continuous site covariate 'habitat', while detection probability can be made dependent an observational covariate 'wind'. Both intercept and slope of the two log- or logistic regressions (for occupancy or expected abundance, respectively, and for detection) are simulated as draws from a normal distribution with  mean and standard deviation that can be selected using function arguments.

Specifically, the data are simulated under the following linear models:

\bold{(1) for \code{type = "det/nondet"} (i.e., community occupancy)}
\tabular{ll}{
(occupancy (psi) and detection (p) for site i, replicate j and species k) \tab \cr
\code{psi[i,k] <- plogis(beta0[k] + beta1[k] * habitat[i]} \tab Occupancy \cr
\code{p[i,j,k] <- plogis(alpha0[k] + alpha1[k] * wind[i,j]} \tab Detection }

\bold{(2) for \code{type = "counts"} (i.e., community count)}
\tabular{ll}{
(exp. abundance (lambda) and detection (p) for site i, rep. j and species k) \tab \cr
\code{lambda[i,k] <- exp(beta0[k] + beta1[k] * habitat[i]} \tab E(N) \cr
\code{p[i,j,k] <- plogis(alpha0[k] + alpha1[k] * wind[i,j]} \tab Detection }

Species-specific heterogeneity in intercepts and slopes is modelled by up to four independent normal distributions (note: no correlation between the intercepts as in Dorazio et al. (2006) or Kéry & Royle (2008))

\bold{(1) for \code{type = "det/nondet"} (i.e., community occupancy)}
\tabular{ll}{
\code{beta0 ~ dnorm(logit(mean.psi), sig.lpsi)} \tab Mean and SD of normal distribution \cr
\code{beta1 ~ dnorm(mu.beta.lpsi, sig.beta.lpsi)} \tab \cr
\code{alpha0 ~ dnorm(logit(mean.p), sig.lp)} \tab \cr
\code{alpha1 ~ dnorm(mu.beta.lp, sig.beta.lp)} \tab }

\bold{(2) for \code{type = "counts"} (i.e., community count)}
\tabular{ll}{
\code{beta0 ~ dnorm(log(mean.lambda), sig.loglam)} \tab Mean and SD of normal distribution \cr
\code{beta1 ~ dnorm(mu.beta.loglam, sig.beta.loglam)} \tab \cr
\code{alpha0 ~ dnorm(logit(mean.p), sig.lp)} \tab \cr
\code{alpha1 ~ dnorm(mu.beta.lp, sig.beta.lp)} \tab }

}
\value{
A list with the arguments supplied and the following additional elements:

\bold{(1) for \code{type = "det/nondet"} (i.e., community occupancy)}

  \item{psi}{\code{psi[i, k]} occupancy probability for site i and species k}
  \item{p}{\code{p[i, j, k]} detection probability for site i, occasion j and species k}
  \item{z}{\code{z[i, k]} true presence/absence for site i and species k}
  \item{z.obs}{\code{z.obs[i, k]} observed presence/absence for site i and species k}
  \item{y.all}{\code{y.all[i, j, k]} detection/nondetection data for site i, occasion j and species k for all species}
  \item{y.obs}{\code{y.obs[i, j, k]} detection/nondetection data for site i, occasion j and species k for observed species}
  \item{y.sum.all}{detection frequency for all species}
  \item{y.sum.obs}{detection frequency for observed species}
  \item{Ntotal.fs}{finite sample (or conditional) species richness}
  \item{Ntotal.obs}{observed species richness}
  \item{S.true}{true number of species occuring at each site}
  \item{S.obs}{observed number of species occuring at each site}
  
\bold{(2) for \code{type = "counts"} (i.e., community count)}
 
  \item{lambda}{\code{lambda[i, k]} expected number of individuals for site i and species k}
  \item{p}{\code{p[i, j, k]} detection probability for site i, occasion j and species k}
  \item{N}{\code{N[i, k]} true number present for site i and species k}
  \item{y.all}{\code{y.all[i, j, k]} number observed for site i, occasion j and species k for all species}
  \item{y.obs}{\code{y.obs[i, j, k]} number observed for site i, occasion j and species k for observed species}
  \item{ymax.obs}{\code{ymax.obs[i, k]} maximum nuumber observed for site i and species k}
  \item{Ntotal.fs}{finite sample (or conditional) species richness}
  \item{Ntotal.obs}{observed species richness}
   
}
\references{
Dorazio et al (2006) Estimating species richness and accumulation by modeling species occurrence and detectability. \emph{Ecology} 87, 842-854

Kéry & Royle (2008) Hierarchical Bayes estimation of species richness and occupancy in spatialy replicated surveys. \emph{J. Appl. Ecol.} 45, 589-598

Royle & Dorazio (2005) 

Yamaura et al. (2012) Biodiversity of man-made open habitats in an underused country: a class of multispecies abundance models for count data. \emph{Biodivers. Conserv.} 21, 1365-1380.

Kéry & Royle (2016) \emph{Applied Hierarachical Modeling in Ecology} Section 11.2

}
\author{
Community occupancy model code partly based on code by Richard Chandler.

}

\examples{
##
}
