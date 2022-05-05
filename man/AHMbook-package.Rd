\name{AHMbook-package}
\alias{AHMbook-package}
\alias{AHMbook}
\docType{package}
\encoding{UTF-8}
\title{
Functions and data for the Book \dQuote{Applied Hierarchical Modeling in Ecology} Volumes 1 and 2
}
\description{
Provides functions to simulate data sets from hierarchical ecological models, including all the simulations described by Marc Kéry and Andy Royle in the two-volume publication, \emph{Applied Hierarchical Modeling in Ecology: analysis of distribution, abundance and species richness in R and BUGS}, Academic Press (Vol 1, 2016; Vol 2, 2021), plus new models developed after publication of the books.

It also has all the utility functions and data sets needed to replicate the analyses shown in the books.
}

\section{SIMULATION FUNCTIONS}{

  \subsection{Introduction}{ \describe{
    \item{\code{\link{sim.fn}}}{ Simulate a homogeneous Poisson point process and illustrate the fundamental relationships between intensity, abundance and occurrence (AHM1 - section 1.1)}
    \item{\code{\link{data.fn}}}{ Simulate count data that are replicated in space and in time according to the binomial N-mixture model of Royle (2004) (this is for much simpler cases than is possible with function \code{simNmix} in Chapter 6 below) (AHM1 - 4.3)}
  }  }

  \subsection{Abundance from Counts: Binomial N-Mixture models}{ \describe{
    \item{\code{\link{simNmix}}}{ Simulate count data and individual detection histories for binomial and multinomial mixture models respectively under a wide range of conditions (AHM1 - 6.9.3)}
    \item{\code{\link{simpleNmix}}}{ Simulate count data under a very simple version of the binomial mixture model, with time for space substitution (AHM1 - 6.12)}
    \item{\code{\link{playRN}}}{ Play Royle-Nichols (RN) model: generate count data under the binomial N-mixture model of Royle (2004),  then 'degrade' the data to detection/nondetection and fit the RN model using \pkg{unmarked} and estimate site-specific abundance (AHM1 - 6.13.1)   }
  }  }

  \subsection{Abundance from Hierarchical Distance Sampling}{ \describe{
    \item{\code{\link{sim.ldata}}}{ Simulate data under a non-hierarchical line transect distance sampling model (AHM1 - 8.2.3)}
    \item{\code{\link{sim.pdata}}}{ Simulate data under a non-hierarchical point transect (= point count) distance sampling model (AHM1 - 8.2.5.1)}
    \item{\code{\link{simHDS}}}{ Simulate data under a hierarchical distance sampling protocol (line or point) (AHM1 - 8.5.1)}
    \item{\code{\link{simHDSpoint}}}{ A simplified version of \code{simHDS} for point transects only.}
    \item{\code{\link{simHDSg}}}{ Simulate data under a hierarchical distance sampling (HDS) protocol with groups (AHM1 - 9.2.1)}
    \item{\code{\link{simHDStr}}}{ Simulate data under a time-removal/distance sampling design (AHM1 - 9.3.2)}
    \item{\code{\link{simHDSopen}}}{ Simulate open hierarchical distance sampling data (AHM1 - 9.5.4)}
    \item{\code{\link{issj.sim}}}{ Simulate data under the open distance sampling protocol for the Island Scrub Jays (AHM1 - 9.7.1)}
    \item{\code{\link{sim.spatialDS}}}{ Simulate data under a basic spatial distance sampling model  (AHM1 - 9.8.3)}
    \item{\code{\link{sim.spatialHDS}}}{ Simulate data under a spatial hierarchical distance sampling model  (AHM1 - 9.8.5)}
    \item{\code{\link{simIDS}}}{ Simulate data for an integrated distance sampling, point count and occupancy study}
  }  }

  \subsection{Static Occurrence using Site-Occupancy Models}{ \describe{
    \item{\code{\link{simOcc}}}{ Simulate detection/nondetection data under static occupancy models under a wide range of conditions (AHM1 - 10.5)}
    \item{\code{\link{simOccCat}}}{ As above, but allows simulation of categorical covariates}
    \item{\code{\link{sim3Occ}}}{ Simulate detection/nondetection data under a static 3-level occupancy model (AHM1 - 10.10)}
    \item{\code{\link{simOccttd}}}{ Simulate 'timing data' under a static time-to-detection occupancy design (AHM1 - 10.12.1)}
    \item{\code{\link{wigglyOcc}}}{ Simulate detection/nondetection data under a static occupancy model with really wiggly covariate relationships in occupancy and detection probability (AHM1 - 10.14)}
  } }

  \subsection{Hierarchical Models for Communities}{ \describe{
    \item{\code{\link{simComm}}}{ Simulate detection/nondetection or count data under a community occupancy or abundance model respectively (AHM1 - 11.2)}
  } }

  \subsection{Relative Abundance Models for Population Dynamics}{ \describe{
    \item{\code{\link{simNpC}}}{ Simulate data on abundance (N), detection probability (p) and resulting counts (C) under a counting process with imperfect detection  (AHM2 - 1.2)}
    \item{\code{\link{simPOP}}}{ Simulate count data under a demographic state-space, or Dail-Madsen, model (no robust design) (AHM2 - 1.7.1)}
    \item{\code{\link{simPH}}}{ Simulate count data with phenological curves within a year (AHM2 - 1.8.1)}
  } }

  \subsection{Modeling Population Dynamics with Count Data}{ \describe{
    \item{\code{\link{simDM0}}}{ Simulate count data from a Dail-Madsen model under a robust design, no covariates (AHM2 - 2.5.1)}
    \item{\code{\link{simDM}}}{  Simulate count data from a Dail-Madsen model under a robust design, with covariates (AHM2 - 2.5.5)}
    \item{\code{\link{simMultMix}}}{  Simulate \dQuote{removal} count data from a multinomial-mixture model (AHM2 - 2.7.1)}
    \item{\code{\link{simFrogDisease}}}{  Simulate detection data for diseased frogs (AHM2 - 2.9.1)}
  } }

  \subsection{Hierarchical Models of Survival}{ \describe{
    \item{\code{\link{simCJS}}}{ Simulate individual capture history data under a Cormack-Jolly-Seber (CJS) survival model (AHM2 - 3.2.2)}
  } }

  \subsection{Dynamic Occupancy Models}{ \describe{
    \item{\code{\link{simDynocc}}}{ Simulate detection/nondetection data under a dynamic occupancy model under a wide range of conditions (AHM2 - 4.4)}
    \item{\code{\link{simDemoDynocc}}}{ Simulate detection/nondetection data under a demographic dynamic occupancy model (AHM2 - 4.12)}
  } }

  \subsection{Dynamic Community Models}{ \describe{
    \item{\code{\link{simDCM}}}{ Simulate detection/nondetection data under a general dynamic community model (site-occupancy variant)  (AHM2 - 5.2)}
  } }

  \subsection{Spatial Models of Distribution and Abundance}{ \describe{
    \item{\code{\link{simDynoccSpatial}}}{ Simulate detection/nondetection data under a dynamic occupancy model with spatial covariate and spatial autocorrelation (AHM2 - 9.6.1.1). See also \code{\link{simDynoccSpatialData}}}
    \item{\code{\link{simExpCorrRF}}}{ Simulate data from a Gaussian random field with negative exponential correlation function (AHM2 - 9.2)}
    \item{\code{\link{simOccSpatial}}}{ Simulate detection/nondetection data under a spatial, static occupancy model for a real landscape in the Bernese Oberland (Swiss Alps) (AHM2 - 9.2)}
    \item{\code{\link{simNmixSpatial}}}{ Simulate counts under a spatial, static binomial N-mixture model for a real landscape in the Bernese Oberland (Swiss Alps) (AHM2 - 9.2)}
  } }

  \subsection{Integrated Models for Multiple Types of Data}{ \describe{
    \item{\code{\link{simPPe}}}{ Simulate a spatial point pattern in a heterogeneous landscape and show aggregation to abundance and occurrence ('e' for educational version)  (AHM2 - 10.2)}
    \item{\code{\link{simDataDK}}}{ Simulate data for an integrated species distribution model (SDM) of Dorazio-Koshkina (AHM2 - 10.6.1)}
  } }

  \subsection{Spatially Explicit Distance Sampling}{ \describe{
    \item{\code{\link{simSpatialDSline}}}{ Simulate line transect distance sampling data with spatial variation in density  (AHM2 - 11.5)}
    \item{\code{\link{simSpatialDSte}}}{ Simulate data for replicate line transect distance sampling surveys with spatial variation in density and temporary emigration (AHM2 - 11.8.1)}
    \item{\code{\link{simDSM}}}{ Simulate line transect data for density surface modeling (AHM2 - 11.10.1)}
  } }
}

\section{DATA SETS}{
  \describe{
  \item{\code{\link{BerneseOberland}}}{ Landscape data for the Bernese Oberland around Interlaken, Switzerland (AHM2 - 9.2)}
  \item{\code{\link{crestedTit}}}{ Crested Tit data from the Swiss Breeding Bird Survey MHB (Monitoring Häufige Brutvögel) for 1999 to 2015 (AHM2 - 1.3)}
  \item{\code{\link{cswa}}}{ Chestnut-sided Warbler data for point counts and spot-mapping from White Mountain National Forest (AHM2 - 2.4.3)}
  \item{\code{\link{crossbillAHM}}}{ Crossbill data from the Swiss Breeding Bird Survey for 2001 to 2012 (AHM2 - 4.9)}
  \item{\code{\link{dragonflies}}}{ Toy data set used in AHM1 - 3.1}
  \item{\code{\link{duskySalamanders}}}{ Counts of juvenile vs adult salamanders over 7 years (AHM2 - 2.9.2)}
  \item{\code{\link{EurasianLynx}}}{ Data for Eurasian Lynx in Italy and Switzerland (AHM2 - 7.3.2)}
  \item{\code{\link{Finnmark}}}{ Data from surveys of birds in Finnmark in NE Norway (AHM2 - 5.7)}
  \item{\code{\link{FrenchPeregrines}}}{ Detection data for peregrines in the French Jura (AHM2 - 4.11)}
  \item{\code{\link{greenWoodpecker}}}{ Count data for Green Woodpeckers in Switzerland from the MHB (AHM2 - 2.2)}
  \item{\code{\link{HubbardBrook}}}{ Point count data for warblers from Hubbard Brook, New Hampshire (AHM2 - 8.2)}
  \item{\code{\link[unmarked]{jay}}}{ The European Jay data set (from the MHB) is now included in \pkg{unmarked} (AHM1 - 7.9)}
  \item{\code{\link{MesoCarnivores}}}{ Camera trap data for 3 species of meso-carnivores (AHM2 - 8.2)}
  \item{\code{\link{MHB2014}}}{ Complete data from the Swiss Breeding Bird Survey MHB (Monitoring Häufige Brutvögel) for the year 2014 (AHM1 - 11.3)}
  \item{\code{\link{spottedWoodpecker}}}{ Data for Middle Spotted Woodpeckers in Switzerland (AHM2 - 4.11.2)}
  \item{\code{\link{SwissAtlasHa}}}{ A 1ha-scale subset of the count data from the Swiss Breeding Bird Atlas (AHM2 - 8.4.2)}
  \item{\code{\link{SwissEagleOwls}}}{ Territory-level, multi-state detection/nondetection data for Eagle Owls in Switzerland (AHM2 - 6.4)}
  \item{\code{\link{SwissMarbledWhite}}}{ Data from the Biodiversity Monitoring Program (LANAG) in the Swiss Canton of Aargau for Marbled White butterfly (AHM2 - 1.8.2)}
  \item{\code{\link{SwissSquirrels}}}{ Count data for Red Squirrels in Switzerland from the Swiss breeding bird survey MHB (AHM1 - 10.9)}
  \item{\code{\link{SwissTits}}}{ Data for 6 species of tits in Switzerland from from the Swiss breeding bird survey MHB  during 2004 to 2013 (AHM1 - 6.13.1)}
  \item{\code{\link{treeSparrow}}}{ Data for Tree Sparrows in Alaska (AHM2 - 11.8.4)}
  \item{\code{\link{ttdPeregrine}}}{ Time-to-detection data for Peregrines (AHM1 - 10.12.2)}
  \item{\code{\link{UKmarbledWhite}}}{ Data from the UK Butterfly Monitoring Scheme (UKBMS) for Marbled White butterfly (AHM2 - 1.8.2)}
  \item{\code{\link{wagtail}}}{ Distance sampling data for Yellow Wagtails in The Netherlands (AHM1 - 9.5.3)}
  \item{\code{\link{waterVoles}}}{ Detection/nondetection data for the Mighty Water Vole of Scotland (AHM2 - 7.2.2)}
  \item{\code{\link{wigglyLine}}}{ Coordinates for a wiggly transect line (AHM2 - 11.9)}
  \item{\code{\link{willowWarbler}}}{ Capture-history (survival) data for Willow Warblers in Britain (AHM2 - 3.4.1)}
} }

\section{UTILITY FUNCTIONS}{
  \describe{
  \item{\code{\link{ppc.plot}}}{ Plot results from posterior predictive checks in section AHM1 - 6.8, for a fitted binomial N-mixture model object with JAGS}
  \item{\code{\link{plot_Nmix_resi}}}{ Do diagnostic plots for one binomial N-mixture model fitted with all three mixture distributions currently available in unmarked: Poisson, negative binomial and zero-inflated Poisson (AHM1 - 6.9.3)}
  \item{\code{\link{map.Nmix.resi}}}{ Produce a map of the residuals from a binomial N-mixture model (see Section AHM1 - 6.9.3)}
  \item{\code{\link{instRemPiFun}}, \code{\link{crPiFun}}, \code{\link{crPiFun.Mb}}, \code{\link{MhPiFun}}}{ Define the relationship between the multinomial cell probabilities and the underlying detection probability parameters (i.e., a pi function) in various designs (AHM1 - 7.8 and AHM2 - Chapter 2)}
  \item{\code{\link{spline.prep}}}{ Prepare input for BUGS model when fitting a spline for a covariate (AHM1 - 10.14)}
  \item{\code{\link{graphSSM}}}{ Plot trajectories of counts and latent abundance from a fitted Gaussian state-space model (AHM2 - 1.6.1)}
  \item{\code{\link{ch2marray}}}{ Convert capture history data to the m-array aggregation (AHM2 - 3.4.1)}
  \item{\code{\link{valid_data}}}{ Partial validation of simulated data with false positives (AHM2 - 7.6.2)}
  \item{\code{\link{getLVcorrMat}}}{ Compute the correlation matrix from an analysis of a latent variable occupancy or binomial N-mixture model (AHM2 - 8.4.2)}
  \item{\code{\link{zinit}}}{ Generate starting values for fitting survival models (introduced in AHM2 - 3.2.3).}
  \item{\code{\link{standardize}}}{ Standardize covariates to mean 0, SD 1.}
  \item{\code{\link{fitstats}}, \code{\link{fitstats2}}}{ Calculate fit-statistics used in parboot GOF tests throughout the book (eg, Sections AHM1 - 7.5.4, AHM1 - 7.9.3, AHM2 - 2.3.3)}
  \item{\code{\link{e2dist}}}{ Compute a matrix of Euclidean distances}
  \item{\code{\link{image_scale}}}{ Draw scale for image (introduced in chapter AHM1 - 9.8.3) }
  \item{\code{\link{bigCrossCorr}}}{ Report cross-correlations above a given threshold }
  \item{\code{\link{Color_Ramps}}}{ Color ramps for use with image or raster plots }
} }


\author{
Marc Kéry, Andy Royle, Mike Meredith
}

\keyword{package}
