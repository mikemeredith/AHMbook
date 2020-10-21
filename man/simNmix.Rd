\name{simNmix}
\alias{simNmix}
\encoding{UTF-8}
\title{
Simulate data for binomial and multinomial mixture models
}
\description{
This very general function generates single-season count data under variants of the binomial N-mixture model of Royle (2004) and of the multinomial N-mixture model of Royle et al (2007).
}
\usage{
simNmix(nsites = 267, nvisits = 3, mean.theta = 1, mean.lam = 2, mean.p = 0.6,
  area = FALSE, beta1.theta = 0, beta2.theta = 0, beta3.theta = 0,
  beta2.lam = 0, beta3.lam = 0, beta4.lam = 0, beta3.p = 0, beta5.p = 0,
  beta6.p = 0, beta.p.survey = 0, beta.p.N = 0, sigma.lam = 0, dispersion = 10,
  sigma.p.site = 0, sigma.p.visit = 0, sigma.p.survey = 0, sigma.p.ind = 0,
  Neg.Bin = FALSE, open.N = FALSE, show.plots = TRUE, verbose = TRUE)
}
\arguments{
  \item{nsites}{
number of sites
}
  \item{nvisits}{
number of visits per site
}
  \item{mean.theta}{
 proportion of sites that can have non-zero abundance in principle: suitability model for zero-inflation
}
  \item{mean.lam}{
Expected abundance at the average value of all abundance covariates (and ignoring random site effects): abundance model
}
  \item{mean.p}{
Expected detection at the average value of all detection covariates (and ignoring all random effects): detection model
}
  \item{area}{
determines area of sites (A), defaults to A=1 (i.e., all identical), but you can supply a vector of site areas of length nsites instead.
}
  \item{beta1.theta}{
coefficient of site covariate 1 in suitability model
}
  \item{beta2.theta}{
coefficient of site covariate 2 in suitability model
}
  \item{beta3.theta}{
coefficient of site covariate 3 in suitability model
}
  \item{beta2.lam}{
coefficient of site covariate 2 in abundance model
}
  \item{beta3.lam}{
coefficient of site covariate 3 in abundance model
}
  \item{beta4.lam}{
coefficient of site covariate 4 in abundance model
}
  \item{beta3.p}{
coefficient of site covariate 3 in detection model
}
  \item{beta5.p}{
coefficient of site covariate 5 in detection model
}
  \item{beta6.p}{
coefficient of site covariate 6 in detection model
}
  \item{beta.p.survey}{
coefficient of survey ('observational') covariate on p
}
  \item{beta.p.N}{
coefficient of centered local population size (log(N+1)) in detection model (i.e., coef. for density-dependent detection prob.)
}
  \item{sigma.lam}{
"Overdispersion SD" in linear predictor of abundance
}
  \item{dispersion}{
 'size' or extra-Poisson dispersion of Negative binomial
}
  \item{sigma.p.site}{
"Overdispersion SD" in linear predictor of detection coming from random site effects
}
  \item{sigma.p.visit}{
 "Overdispersion SD" in linear predictor of detection coming from random visit effects
}
  \item{sigma.p.survey}{
"Overdispersion SD" in linear predictor of detection coming from random site-by-survey effects
}
  \item{sigma.p.ind}{
"Overdispersion SD" in linear predictor of detection coming from random site-by-individual effects
}
  \item{Neg.Bin}{
if FALSE, any overdispersion in abundance is modeled by a Poisson log-normal; if TRUE, abundance overdispersion is modeled by adoption of a Negative binomial distribution for latent N
}
  \item{open.N}{
if TRUE, data are simulated under one specific form of an open population, where N in the first occasion is drawn from the specified mixture distribution and for all further occasions j, we have  N_ij ~ Poisson(N_i(j-1)). With open.N = TRUE, we must have sigma.p.ind = 0, show.plot = FALSE and nvisits >1.
}
  \item{show.plots}{
 if TRUE, plots of the data will be displayed; set to FALSE if you are running simulations.
}
  \item{verbose}{
if TRUE, output will be written to the console.
}
}
\details{
Data are simulated at the level of each individual and individual-specific detection heterogeneity can be included. As a side-effect, individual-specific detection histories are generated and hence, data are also be simulated under the corresponding multinomial N-mixture model.

Broadly, the function can generate data under this most general model:

'Suitability' (zero-inflation) ~ cov1 + cov2 + cov3

Abundance ~ offset + cov2 + cov3 + cov4 + overdispersion

Detection ~ cov3 + cov5 + cov6 + survey.covariate + log(N+1) + eps.site + eps.visit + eps.survey + eps.individual

Overdispersion in abundance is modeled either as a Poisson-log-normal with a normal random site effect in lambda or with a Negative binomial with mean lambda and a 'size', or dispersion, parameter. Variable site areas can be specified to affect abundance as in an offset.

Abundance can be zero-inflated (this is the 'suitability' model). Note that the zero-inflation parameter is called theta here (in unmarked it is called psi). mean.phi is the probability that a site is suitable (i.e., 1 minus the expected proportion of sites with structural zero abundance.

Site covariate 2 can affect both suitability and abundance, while covariate 3 may affect all three levels. Hence, the function permits to simulate the case where a single site covariate affects different levels in the process (e.g., abundance and detection) in opposing directions (as for instance in Kéry, 2008)

Density-dependent detection can be modeled as a logistic-linear effect of local abundance (centered and log(x+1) transformed). Overdispersion in detection is modeled via normal random effects (the eps terms above) specific to sites, visits, surveys or individuals.

Effects of covariates and random-effects factors are modeled as additive on the link scale (log for abundance and logit for suitability and detection).

Data may be generated under one specific open-population model when argument 'open.N' is set to TRUE.
}

\value{
A list with the arguments input and the following additional elements:
 \item{nobs }{The total number of visits}
 \item{site.cov }{An nsites x 6 matrix of values for 6 site covariates}
 \item{survey.cov }{An nsites x nvisits matrix of values for a survey covariate}
 \item{log.lam }{Linear predictor of PLN abundance model including random effects, a vector of length nsites}
 \item{s }{Site suitability indicator, a vector of length nsites}
 \item{N }{Number of individuals at each site, a vector of length nsites}
 \item{p }{Probability of detection, an array with dimensions sites occupied x visits x max(N)}
 \item{DH }{Detection history (1/0), an array with dimensions sites occupied x visits x max(N)}
 \item{N.open }{Number of individuals at each site for open model, a sites x visits matrix}
 \item{C }{Summary of DH: number of individuals detected for each site and visit}
 \item{eta.lam }{Random site effects in lambda, a vector of length nsites, zero if Neg.Bin == TRUE}
 \item{eta.p.site }{Random site effects in p, a vector of length nsites}
 \item{eta.p.visit }{Random visit effects in p, a vector of length nvisits}
 \item{eta.p.survey }{Random survey (= site-by-survey) effects in p, a nsites x nvisits matrix}
 \item{eta.p.ind }{Random individual (= site-by-ind) effects in p (NOT site-ind-visit !), a nsites x max(N) matrix}
 \item{odcN }{Naive overdispersion measure (var/mean) for true abundance (N)}
 \item{odcC }{Naive overdispersion measure (var/mean) for observed counts (C)}
 \item{Ntotal }{Total abundance summed over all sites}
 \item{summax }{The sum of maximum counts over all sites}

}
\references{
Royle, J.A. (2004) N-mixture models for estimating population size from spatially replicated counts. \emph{Biometrics}, 60, 108-115.

Royle, J.A., et al (2007) Hierarchical spatial models of abundance and occurrence from imperfect survey data. \emph{Ecological Monographs}, 77, 465-481.

Kéry, M. (2008) Estimating abundance from bird counts: binomial mixture models uncover complex covariate relationships, \emph{Auk}, 125, 336-345

Kéry, M. & Royle, J.A. (2016) \emph{Applied Hierarchical Modeling in Ecology} AHM1 - 6.5.
}
\author{
Marc Kéry & Andy Royle
}
\examples{
# Generate data with the default arguments and look at the structure:
tmp <- simNmix()
str(tmp)

\donttest{
str(data <- simNmix())                  # Null data-generating model
str(data <- simNmix(mean.theta = 0.60)) # ZIP with 40% structural zeroes
str(data <- simNmix(sigma.lam = 1))     # Poisson-lognormal (PLN) mixture
str(data <- simNmix(Neg.Bin = TRUE))    # Negative-binomial mixture
str(data <- simNmix(mean.theta = 0.6, sigma.lam = 1))  # Zero-inflated PLN
str(data <- simNmix(mean.theta = 0.6, Neg.Bin = TRUE)) # Zero-infl. NegBin
str(data <- simNmix(mean.p = 1))        # Perfect detection (p = 1)
str(data <- simNmix(mean.theta = 0.6, mean.p = 1))     # ZIP with p = 1
str(data <- simNmix(sigma.lam = 1, mean.p = 1))        # PLN with p = 1
}
}
