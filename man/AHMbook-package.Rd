\name{AHMbook-package}
\alias{AHMbook-package}
\alias{AHMbook}
\docType{package}
\encoding{UTF-8}
\title{
Functions and data for the Book 'Applied Hierarchical Modeling in Ecology'
}
\description{
Provides functions and data sets needed to run the code given in the two-volume publication, \emph{Applied Hierarchical Modeling in Ecology: analysis of distribution, abundance and species richness in R and BUGS} by Marc Kéry and Andy Royle, Academic Press (Vol 1, 2016; Vol 2, in preparation).
}
\details{


The functions are listed by chapter below.
}

\section{Chapter 1}{
  \describe{
  \item{\code{\link{sim.fn}}}{ Simulate a Poisson point process and illustrate the relationships between intensity, abundance and occurrence (section 1.1)}
  } }

\section{Chapter 4}{
  \describe{
  \item{\code{\link{data.fn}}}{ Simulate count data that are replicated in space and in time according to the binomial mixture model of Royle (2004) (this is for much simpler cases than is possible with \code{simNmix} below) (4.3)}
  }  }

\section{Chapter 6}{
  \describe{
  \item{\code{\link{ppc.plot}}}{ Plot results from posterior predictive checks in section 6.8, for a fitted model object with JAGS}
  \item{\code{\link{simNmix}}}{ Simulate data for binomial and multinomial mixture models under a wide range of conditions (6.9.3)}
  \item{\code{\link{plot_Nmix_resi}}}{ Do diagnostic plots for one Nmix model fitted with all three mixture distributions currently available in unmarked: Poisson, negative binomial and zero-inflated Poisson (6.9.3)}
  \item{\code{\link{map.Nmix.resi}}}{ Produce a map of the residuals from an Nmix model (see Section 6.9.3)}
  \item{\code{\link{simpleNmix}}}{ Simulate data under a simple version of the binomial mixture model, with space for time substitution (6.12)}
  \item{\code{\link{playRN}}}{ Play Royle-Nichols model: generate replicated count data under the Nmix model of Royle (2004),  then 'degrade' the data to detection/nondetection and fit the RN model using unmarked and estimate site-specific abundance (6.13.1)   }
  } }

\section{Chapter 7}{
  \describe{
  \item{\code{\link{instRemPiFun}}, \code{\link{crPiFun}}, \code{\link{crPiFun.Mb}}, \code{\link{MhPiFun}}}{ Define the relationship between the multinomial cell probabilities and the underlying detection probability parameters (i.e., a pi function) in various designs (7.8 and Chapter 13)}
  } }

\section{Chapter 8}{
  \describe{
  \item{\code{\link{sim.ldata}}}{ Simulate data under a non-hierarchical line transect distance sampling model (8.2.3)}
  \item{\code{\link{sim.pdata}}}{ Simulate data under a non-hierarchical point transect (= point count) distance sampling model (8.2.5.1)}
  \item{\code{\link{simHDS}}}{ Simulate data under a hierarchical distance sampling protocol (line or point) (8.5.1)}
  }  }

\section{Chapter 9}{
  \describe{
  \item{\code{\link{simHDSg}}}{ Simulate data under a hierarchical distance sampling (HDS) protocol with groups (9.2.1)}
  \item{\code{\link{simHDStr}}}{ Simulate data under a time-removal/distance sampling design (9.3.2)}
  \item{\code{\link{simHDSopen}}}{ Simulate open hierarchical distance sampling data (9.5.4)}
  \item{\code{\link{issj.sim}}}{ Simulate data under the open distance sampling protocol for the Island Scrub Jays (9.7.1)}
  \item{\code{\link{sim.spatialDS}}}{ Simulate data under a basic spatial distance sampling model  (9.8.3)}
  \item{\code{\link{sim.spatialHDS}}}{ Simulate data under a spatial hierarchical distance sampling model  (9.8.5)}
  } }

\section{Chapter 10}{
  \describe{
  \item{\code{\link{simOcc}}}{ Simulate data under static occupancy models under a wide range of conditions (10.5)}
  \item{\code{\link{sim3Occ}}}{ Simulate data under a static 3-level occupancy model (10.10)}
  \item{\code{\link{simOccttd}}}{ Simulate data under a static time-to-detection occupancy design (10.12.1)}
  \item{\code{\link{wigglyOcc}}}{ Simulate data under a static occupancy model with really wiggly covariate relationships in occupancy and detection probability (10.14)}
  \item{\code{\link{spline.prep}}}{ Prepare input for BUGS model when fitting a spline for a covariate (10.14)}
  } }

\section{Chapter 11}{
  \describe{
  \item{\code{\link{simComm}}}{ Simulate community occupancy or abundance data (11.2)}
  } }

\section{Chapter 12}{
  \describe{
  \item{\code{\link{simNpC}}}{ Simulate N, p and C (12.2)}
  \item{\code{\link{simPOP}}}{ Simulate demographic state-space model (12.7.1)}
  \item{\code{\link{simPH}}}{ Simulate phenology data (12.8.1)}
  \item{\code{\link{graphSSM}}}{ Plot state-space model parameters (12.6)}
  } }

\section{Chapter 13}{
  \describe{
  \item{\code{\link{simDM0}}}{ Simulate data from a Dail-Madsen model, no covariates (13.5.1)}
  \item{\code{\link{simDM}}}{  Simulate data from a Dail-Madsen model, with covariates (13.5.5)}
  \item{\code{\link{simMultMix}}}{  Simulate "removal" data from a multinomial-mixture model (13.7.1)}
  } }

\section{Chapter 14}{
  \describe{
  \item{\code{\link{simCJS}}}{ Simulate data for CJS survival models (14.2.2)}
  \item{\code{\link{ch2marray}}}{ Convert capture history data to marray (14.4.1)}
  } }

\section{Chapter 15}{
  \describe{
  \item{\code{\link{simDynocc}}}{ Simulate data under a dynamic occupancy model under a wide range of conditions}
  \item{\code{\link{simDemoDynocc}}}{ Simulate data under a demographic dynamic occupancy model}
  } }

\section{Chapter 16}{
  \describe{
  \item{\code{\link{simDCM}}}{ Simulate detection/nondetection data under a general dynamic community (site-occupancy) model}
  } }

\section{Chapter 18}{
  \describe{
  \item{\code{\link{valid_data}}}{ Partial validation of simulated data with false positives}
  \item{\code{\link{modSelFP}}}{ A stop-gap function for \code{\link[unmarked]{modSel}} for fitted false positive models }
  } }

\section{Chapter 20}{
  \describe{
  \item{\code{\link{spatial.exp}}}{ Generate Gaussian random field with negative exponential correlation}
  \item{\code{\link{simDynoccSpatial}}}{ Simulate data under a dynamic occupancy model with spatial covariate and spatial autocorrelation}
  } }

\section{Chapter 21}{
  \describe{
  \item{\code{\link{simPPe}}}{ Simulate a spatial point pattern in a heterogeneous landscape}
  \item{\code{\link{simDataDK}}}{ Simulate data for an integrated Dorazio-Koshkina analysis}
  } }

\section{Chapter 22}{
  \describe{
  \item{\code{\link{sim.spatialDSline}}}{ Simulate data for a standardized line transect}
  \item{\code{\link{sim.spatialDSte}}}{ Simulate data for replicate line transect surveys with temporary emigration}
  } }

\section{DATA SETS}{
  \describe{
  \item{\code{\link{crestedTit}}}{ Crested tit data from the Swiss Breeding Bird Survey for 1999 to 2015 (12.3)}
  \item{\code{\link{cswa}}}{ Chestnut-sided warbler data for point counts and spot-mapping from White Mountain National Forest (13.4.3)}
  \item{\code{\link{crossbillAHM}}}{ Crossbill data from the Swiss Breeding Bird Survey for 2001 to 2012 (15...)}
  \item{\code{\link{dragonflies}}}{ Toy data set used in section 3.1}
  \item{\code{\link{duskySalamanders}}}{ Counts of juvenile vs adult salamanders over 7 years (13.9.2).}
  \item{\code{\link{EurasianLynx}}}{ Data for lynx in Italy and Switzerland (18.3.2)}
  \item{\code{\link{Finnmark}}}{ Data from surveys of birds in Finnmark in NE Norway (chapter 17)}
  \item{\code{\link{FrenchPeregrines}}}{ Detection data for peregrines in the French Jura (chapter 15)}
  \item{\code{\link{greenWoodpecker}}}{ Count data for green woodpeckers in Switzerland (13.2)}
  \item{\code{\link{HubbardBrook}}}{ Point count data for warblers from Hubbard Brook, New Hampshire (19.2)}
  \item{\code{\link[unmarked]{jay}}}{ The European jay data set is now included in \pkg{unmarked} (7.9)}
  \item{\code{\link{MHB2014}}}{ Data from the Swiss Breeding Bird Survey for 2014 (11.3)}
  \item{\code{\link{spottedWoodpecker}}}{ Data for middle spotted woodpeckers in Switzerland (15.11.2)}
  \item{\code{\link{SwissAtlasHa}}}{ Data from the Swiss Bird Atlas at the scale of 1ha (19..)}
 \item{\code{\link{SwissEagleOwls}}}{ Data for Eagle Owls in Switzerland (17..)}
  \item{\code{\link{SwissMarbledWhite}}}{ Data from the Aargau Biodiversity Monitoring Program (LANAG) for marbled white butterflies (12.8.2)}
  \item{\code{\link{SwissSquirrels}}}{ Data for red squirrels in Switzerland (10.9)}
  \item{\code{\link{SwissTits}}}{ Data for 6 species of tits in Switzerland from 2004 to 2013 (6.13.1)}
  \item{\code{\link{ttdPeregrine}}}{ Time-to-detection data for peregrines (10.12.2)}
 \item{\code{\link{UKmarbledWhite}}}{ Data from the UK Butterfly Monitoring Scheme (UKBMS) for marbled white butterflies (15..)}
  \item{\code{\link{wagtail}}}{ Data for yellow wagtails in Netherlands (9.5.3)}
  \item{\code{\link{waterVoles}}}{ Data for water voles from Scotland (18.2.2)}
  \item{\code{\link{willowWarbler}}}{ Data for willow warblers in Britain (14.4.1)}

  } }

\section{UTILITY FUNCTIONS}{
  \describe{
  \item{\code{\link{zinit}}}{ Starting values for fitting survival models (introduced in chapter 15).}
  \item{\code{\link{standardize}}}{ Standardize covariates to mean 0, SD 1.}
  \item{\code{\link{fitstats}}, \code{\link{fitstats2}}}{ Calculate fit-statistics used in parboot GOF tests throughout book (Sections 7.5.4, 7.9.3, 13.3.3)}
  \item{\code{\link{e2dist}}}{ Compute Euclidean distances}
  \item{\code{\link{image_scale}}}{ Draw scale for image (from SCR book, introduced in chapter 9) }
  \item{\code{\link{bigCrossCorr}}}{ Report cross-correlations above a given threshold }
  \item{\code{\link{Color_Ramps}}}{ Color ramps for use with image or raster plots }
  } }


\author{
Marc Kéry, Andy Royle, Mike Meredith

}

\keyword{package}
