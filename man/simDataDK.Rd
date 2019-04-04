\name{simDataDK}
\alias{simDataDK}
\encoding{UTF-8}

\title{
Simulate data for an integrated Dorazio-Koshkina analysis
}
\description{
The function generates a population in a heterogeneous landscape and simulates data from two different sources: (1) presence-only data, and (2) replicate counts in randomly-selected quadrats. This is the scenario for the integrated model described by Dorazio (2014) and Koshkina et al. (2017).

A Poisson point pattern (PPP) with intensity a function of a covariate X and intercept and coefficient beta is simulated on a discrete (pixel-based) approximation of a continuous landscape.

This PPP is then thinned with a pixel-wise thinning probability controlled by a covariate W and coefficients alpha, and then with a landscape-wise random drop-out process to produce a first data set of presence-only kind.

A second data set is simulated by imagining replicated counts conducted in randomly-selected quadrats within the landscape. Detection of individuals is imperfect, with probability of detection controlled by the covariate W and coefficients gamma.
}
\usage{
simDataDK(npix = 100, alpha = c(-1,-1), beta = c(6,0.5),
  drop.out.prop.pb = 0.7, quadrat.size = 4, gamma = c(0,-1.5),
  nquadrats = 250, nsurveys = 3, show.plot = TRUE)
}

\arguments{
  \item{npix}{
number of pixels along each side of square state space (the 'landscape'); the total number of pixels is then npix^2.
}
  \item{alpha}{
coefficients for the relationship: logit(b) = alpha[1] + alpha[2] * W, where b is the sampling detection bias in the presence-only observations.
}
  \item{beta}{
coefficients for the relationship: log(lambda) = beta[1] + beta[2] * X, where lambda is the intensity of the Poisson point process. For computational reasons, each pixel is limited to one animal; if the values of beta require more animals than pixels in the highest intensity regions, an error will occur. %% Error message: "The 'beta' settings result in intensities that are too high in the most intense region (more animals than pixels)."
}
  \item{drop.out.prop.pb}{
proportion of presence-only points at the end that are discarded.
}

  \item{quadrat.size}{
length of the side of quadrats for conducting replicate counts in \emph{pixel} units; note that npix / quadrat.size must yield an integer.
}
  \item{gamma}{
coefficients for the relationship: logit(p) = gamma[1] + gamma[2] * W, where p is the probability of detecting an individual during the count surveys in the quadrats.
}
  \item{nquadrats}{
the number of quadrats selected for the count survey.
}
  \item{nsurveys}{
the number of replicate counts in each quadrat.
}
  \item{show.plot}{
if TRUE, summary plots are displayed.
}
}

\value{
A list with the values of the input arguments and the following additional elements:
 \item{s.area }{the area of the whole landscape = 4}
 \item{s.loc }{2-column vector with the location of each pixel}
 \item{xcov }{values of the 'X' (intensity) covariate}
 \item{wcov }{values of the 'W' (detection) covariate}
 \item{N.ipp }{true number of individuals in the landscape}
 \item{pixel.id.ipp }{pixel ID for each individual in the population}
 \item{loc.ipp }{coordinates for each individual in the population}
 \item{pTrue.ipp }{probability of detection for each individual for presence-only data}
 \item{pixel.id.det1 }{pixel ID for each individual detected before drop-out}
 \item{N.det1 }{number of detections before drop-out}
 \item{det.ipp1 }{coordinates of each individual detected before drop-out}
 \item{pixel.id.det2 }{pixel ID for each individual detected after drop-out}
 \item{N.det2 }{number of detections after drop-out}
 \item{det.ipp2 }{coordinates of each individual detected after drop-out}
 \item{pcount }{probability of detection during count surveys, varies by quadrat}
 \item{fullCountData }{matrix with rows for each quadrat, columns for ID, x and w coords, true N, and 3 replicate counts}
 \item{countData }{as above, but rows for quadrats sampled only}
 \item{s }{a Raster Stack with layers for 'X' and 'W', lambda and pTrue}
 \item{squad }{a Raster Stack corresponding to the quadrats, with mean 'X' and 'W' and true abundance, 'N'}
}
\references{
Dorazio, R.M. (2014) Accounting for imperfect detection and survey bias in statistical analysis of presence-only data. \emph{Global Ecology and Biogeography}, 23, 1472-1484.

Koshkina, V., Wang, Y., Gordon, A., Dorazio, R.M., White, M., & Stone, L. (2017) Integrated species distribution models: combining presence-background data and site-occupany data with imperfect detection. \emph{Methods in Ecology and Evolution}, 8, 420-430.

Kéry & Royle (2019) \emph{Applied Hierarchical Modeling in Ecology} Chapter 21
}
\author{
Marc Kéry, based on the code written by Dorazio (2014) and adapted by Koshkina et al. (2017).
}

\examples{
# Run the function with default values and look at the output
str(tmp <- simDataDK(), 1)  # use str(., max.level=1) to limit the amount of output.

str(tmp <- simDataDK(show.plot=FALSE), 1)  # no plots

str(tmp <- simDataDK(npix = 500), 1)       # much larger landscape

str(tmp <- simDataDK(alpha = c(-1,1)), 1)  # positive effect of W on bias rate parameter b

str(tmp <- simDataDK(beta = c(6, 0.5)), 1) # lower density

str(tmp <- simDataDK(drop.out.prop = 0), 1)# No final uniform thinning ("drop out")

str(tmp <- simDataDK(beta = c(6, 1)), 1)    # steeper gradient of habitat suitability
}

