\name{AHMbook-package}
\alias{AHMbook-package}
\alias{AHMbook}
\docType{package}
\encoding{UTF-8}
\title{
Functions and data for the Book 'Applied Hierarchical Modeling in Ecology'
}
\description{
Provides functions and data sets needed to run the code given in the two-volume publication, \emph{Applied Hierarchical Modeling in Ecology: analysis of distribution, abundance and species rishness in R and BUGS} by Marc Kéry and Andy Royle, Academic Press (Vol 1, 2016; Vol 2, in preparation).
}
\details{


The functions are listed by chapter below.
}

\section{Chapter 1}{
 
  \tabular{ll}{
  \code{\link{sim.fn}} \tab Simulation of a Poisson point process and illustration of the relationships between the intensity, abundance and occurrence} }
  
\section{Chapter 4}{
 
  \tabular{ll}{
  \code{\link{data.fn}} \tab Simulation of count data that are replicated in space and in time according to the binomial mixture model of Royle (2004) (this is for much simpler cases than what is possible with simNmix, see below)} }
  
\section{Chapter 6}{
 
  \tabular{ll}{
  \code{\link{simNmix}} \tab Simulation of data for binomial and multinomial mixture models under a wide range of conditions \cr
  \code{\link{ppc.plot}} \tab Plot results from posterior predictive check in AHM1, section 6.8, for a fitted model object with JAGS \cr
  \code{\link{plot_Nmix_resi}} \tab Do diagnostic plots for one Nmix model fitted with all three mixture distributions currently available in unmarked: Poisson, negative binomial and zero-inflated Poisson (as per Section 6.9.3 in the AHM1 book) \cr
  \code{\link{map.Nmix.resi}} \tab  Produce a map of the residuals from an Nmix model (see Section 6.9.3 in the AHM1 book) \cr
  \code{\link{simpleNmix}} \tab  Simulation of data under a simple version of the binomial mixture model, with space for time substitution (see AHM1 section 6.12) \cr
  \code{\link{playRN}} \tab Play Royle-Nichols model: generate replicated count data under the Nmix model of Royle (2004),  then 'degrade' the data to detection/nondetection and fit the RN model using unmarked and estimate site-specific abundance (see AHM1 section 6.13.1)   } }
  
\section{Chapter 7}{
 
  \tabular{ll}{
  \code{\link{fitstats}} \tab Returns three fit-statistics (used in parboot GOF tests throughout book) (used, among others, in Chapter 7, e.g., Section 7.5.4)\cr
  \code{\link{fitstats2}} \tab Returns several fit-statistics (used in parboot GOF tests in section 7.9.3)\cr
  \code{\link{instRemPiFun}} \tab Define the relationship between the multinomial cell probabilties and the underlying detection probability parameters (i.e., a pi function) in a removal design with unequal length of sampling periods. \cr
  \code{\link{crPiFun}} \tab Define a pi function for capture-recapture design with 3 surveys (introduced in Section 7.8) \cr
  \code{\link{crPiFun.Mb}} \tab Define a pi function for capture-recapture design with 3 surveys and behavioural response (introduced in Section 7.8.2) \cr
  \code{\link{MhPiFun}} \tab Define a pi function for model with individual detection heterogeneity (introduced in Section 7.8.3) } }
  
\section{Chapter 8}{
 
  \tabular{ll}{
  \code{\link{sim.ldata}} \tab Simulation of data under a non-hierarchical line transect distance sampling model (introduced in AHM1 Section 8.2.3)\cr
  \code{\link{sim.pdata}} \tab Simulation of data under a non-hierarchical point transect (= point count) distance sampling model (introduced in Section 8.2.5.1)\cr
  \code{\link{simHDS}} \tab Simulation of data under a hierarchical distance sampling protocol (line or point) (introduced in AHM1 Section 8.5.1)} }
  
\section{Chapter 9}{
 
  \tabular{ll}{
  \code{\link{simHDSg}} \tab Simulation of data under a hierarchical distance sampling (HDS) protocol with groups (introduced in AHM1 Section 9.2.1)\cr
  \code{\link{simHDStr}} \tab Simulation of data under a time-removal/distance sampling design (introduced in Section 9.3.2)\cr
  \code{\link{simHDSopen}} \tab Simulate open hierarchical distance sampling data \cr
  \code{\link{issj.sim}} \tab Simulation of data under the open distance sampling protocol for the Island Scrub Jays (introduced in Section 9.7.1) \cr
  \code{\link{sim.spatialDS}} \tab Simulation of data under a spatial hierarchical distance sampling model  (introduced in Section 9.8.3)\cr
  \code{\link{pcount.spHDS}} \tab Fits a spatial hierarchical distance sampling model using a hacked version of the pcount function in unmarked (introduced in Section 9.8.4)} }

\section{Chapter 10}{
 
  \tabular{ll}{
  \code{\link{simOcc}} \tab Simulation of data under static occupancy models under a wide range of conditions (introduced in AHM1, Section 10.5)\cr
  \code{\link{sim3Occ}} \tab Simulation of data under a static 3-level occupancy model (introduced in AHM1, section 10.10)\cr
  \code{\link{simOccttd}} \tab Simulation of data under a static time-to-detection occupancy design (Garrard et al., Austral Ecology, 2008, also see Bornand et al., MEE, 2014; introduced in AHM1, section 10.12.1)\cr
  \code{\link{wigglyOcc}} \tab Simulation of data under a static occupancy model with really wiggly covariate relationships in occupancy and detection probability (introduced in AHM1, section 10.14)\cr
  \code{\link{spline.prep}} \tab Prepare input for BUGS model when fitting a spline for a covariate (introduced in AHM1, section 10.14)} }

\section{Chapter 11}{
 
  \tabular{ll}{
  \code{\link{simComm}} \tab Simulation community occupancy or abundance data } }
  
\section{Chapter 15}{
 
  \tabular{ll}{
  \code{\link{simDynocc.fn}} \tab Simulation of data under a dynamic occupancy model under a wide range of conditions} }
  
\section{Chapter 20}{
 
  \tabular{ll}{
  \code{\link{spatial.exp}} \tab ...} }
  
\section{UTILITY FUNCTIONS}{

  \tabular{ll}{
  \code{\link{e2dist}} \tab Compute Euclidean distances \cr
  \code{\link{image.scale}} \tab Draw scale for image (from SCR book, introduced in AHM1 chapter 9) }
}

\author{
Marc Kéry, Andy Royle, Mike Meredith

}

\keyword{ package }
