# Code from Marc's Word file "new simHDS function with availability.docx"
# See emails 2022-02-16 ff.

# 1 A new data simulation function for HDS + PC data
# We now package into a function code written by Ken. All distance units will be thought of as being in metres

# Changes by Mike:
# Change function name to 'simIDS'.
# Use new 'simHDSpoint' function.
# Use density instead of abundance, change mean.lambda to mean.density, beta.lam to beta.density.
# Reorder arguments by data type (DS, PC, OC); use upper case HDS, PC, OC for all arguments.
# Skip PC or OC data generation if nsites_PC or nsites_OC = 0


# ------------ Start of function definition ---------------
simIDS <- function(mean.density = 1, beta.density = 1, mean.phi = 0.14, beta.phi = 0,
    # for distance sampling
    nsites_HDS = 1000,  sigHDS = 100, maxDist_HDS = 200, nbins = 4, range.dur.HDS = c(5, 5),
    # for point counts
    nsites_PC = 10000, sigPC = 70, maxDist_PC = 500, range.dur.PC = c(3, 30),
    # for occupancy = detection/nondetection
    nsites_OC = 5000, sigOC = sigPC, maxDist_OC = maxDist_PC, range.dur.OC = range.dur.PC,
    show.plots = TRUE) {
  #
  # Generates HDS data and PC data for an integrated distance sampling
  # (IDS) model with shared density and availability processes,
  # but possibly different detection/perceptability process
  # (i.e., different detection functions).
  #
  # All lengths are in meters, densities in animals per hectare.
  #
  # Arguments:
  # Ecological model:
  # mean.density: expected density, animals per ha.
  # beta.density: coefficient of expected density on habitat covariate
  # mean.phi: singing or activity rate phi, so that under the model of
  # Solymos et al. (2013) availability probability
  # theta = 1 – exp(- duration * mean.phi).
  # beta.phi: coefficient of log(singing rate) on some covariate (LATER PERHAPS)

  # Distance sampling
  # nsites_HDS: number of sites with HDS protocol
  # sigHDS: scale parameter sigma in the half-normal detection function
  #    at the HDS sites
  # nbins:   Number of distance bins
  # maxDist_HDS: Truncation distance of observations in the HDS protocol;
  #    any observations beyond this are discarded.
  # range.dur.HDS: range of the survey durations in minutes for the HDS
  #    surveys. By default, this is set to c(5, 5) for the HDS data,
  #    to enforce constancy (with 5 min distance sampling)

  # Point counts
  # nsites_PC: number of sites with simple point count (PC) protocol
  # sigPC: scale parameter sigma in the half-normal detection function
  #    at the PC sites
  # maxDist_PC: Maximum distance from the observer at which animals can be detected.
  # range.dur.PC: range of the survey durations in minutes for the PC
  #    surveys. Means will draw from uniform distribution on those bounds.

  # Detection/nondetection
  # nsites_OC: number of sites with Detection/nondetection (occupancy) protocol
  # sigOC: scale parameter sigma in the half-normal detection function
  #    at the OC sites
  # maxDist_OC: Maximum distance from the observer at which animals can be detected.
  # range.dur.OC: range of the survey durations in minutes for the OC
  #    surveys. Means will draw from uniform distribution on those bounds.

  # Written by Ken Kellner and Marc Kéry, February 2022
  # Desecrated by Mike Meredith, on-going


    # Distance sampling data
  # ----------------------
  ## Simulate a regular distance sampling data set
  dat1.raw <- simHDSpoint(nsites = nsites_HDS,
    mean.density = mean.density, beta.density = beta.density, mean.sigma = sigHDS,
    beta.sigma = 0, B = maxDist_HDS, discard0 = FALSE, show.plots = show.plots)
  # str(dat1.raw)

  # Re-format to put into an `unmarkedFrameDS` object:
  # Convert long-form to y matrix
  db <- seq(0, maxDist_HDS, length.out=(nbins+1)) # Distance bins, equal widths

  tmp <- with(dat1.raw, tapply(data[,2], data[,1], function(x) hist(x, breaks=db, plot=FALSE)$counts)) # returns a list
  y_hds <- do.call(rbind, tmp)  # convert list to a matrix

  # Site covariates
  sc_hds <- data.frame(habitat=dat1.raw$habitat)

  # Create unmarked frame for HDS data set:
  #   'Null data set' has perfect availability
  umf_hds0 <- unmarkedFrameDS(y=y_hds, siteCovs=sc_hds, survey="point",
              dist.breaks=db, unitsIn="m")
  # summary(umf_hds0)

  # Add availability
  # Simulate survey durations for HDS data
  dds <- runif(nsites_HDS, min(range.dur.HDS), max(range.dur.HDS)) # HDS

  # Linear model for the singing/activity rate phi
  phi <- mean.phi              # possibly later also add covariate

  ### Simulate availability process on existing datasets
  # (1) Model of Solymos et al. 2013
  # HDS
  pdds <- 1-exp(-1*dds*phi) # Compute availability prob.
  umf_hds1 <- umf_hds0           # Copy unmarked data frame
  tmp <- umf_hds0@y
  umf_hds1@y <- array(rbinom(length(tmp), tmp, pdds), dim=dim(tmp))

  # Point count data
  # ----------------
  if(nsites_PC > 0) {
    ## Simulate Point Count data set (raw, without availability so far)
    dat2.raw <- simHDSpoint(nsites = nsites_PC,
      mean.density = mean.density, beta.density = beta.density, mean.sigma = sigPC,
      beta.sigma = 0, B = maxDist_PC, discard0 = FALSE, show.plots = show.plots)

    # Re-format to put in a unmarkedFramePCount object
    y_pc <- matrix(dat2.raw$counts, ncol=1)

    # Site covariates
    sc_pc <- data.frame(habitat=dat2.raw$habitat)

    # Create unmarked frame for PC data set:
    #   'Null data set' with perfect availability
    umf_pc0 <- unmarkedFramePCount(y=y_pc, siteCovs=sc_pc)
    # summary(umf_pc0)

    # Add availability
    # Simulate survey durations for PC data
    dpc <- runif(nsites_PC, min(range.dur.PC), max(range.dur.PC))

    # Linear model for the singing/activity rate phi
    phi <- mean.phi              # possibly later also add covariate

    # Simulate availability process on existing datasets
    pdpc <- 1-exp(-1*dpc*phi) # Compute availability prob.
    umf_pc1 <- umf_pc0             # Copy unmarked data frame
    # tmp <- umf_pc0@y
    # umf_pc1@y <- array(rbinom(length(tmp), tmp, pdpc), dim=dim(tmp))
    umf_pc1@y <- matrix(rbinom(nsites_PC, umf_pc0@y, pdpc), ncol=1)
  }

  # Detection/nondetection data
  # ---------------------------
  if(nsites_OC > 0) {
    ## Simulate occupancy data set (raw, without availability so far)
    dat3.raw <- simHDSpoint(nsites = nsites_OC,
      mean.density = mean.density, beta.density = beta.density, mean.sigma = sigOC,
      beta.sigma = 0, B = maxDist_OC, discard0 = FALSE, show.plots = show.plots)

    # Re-format to put in a unmarkedFramePCount object
    tmp <- dat3.raw$counts > 0
    y_oc <- matrix(as.numeric(tmp), ncol=1)

    # Site covariates
    sc_oc <- data.frame(habitat=dat3.raw$habitat)

    # Create unmarked frame for OC data set:
    #   'Null data set' with perfect availability
    umf_oc0 <- unmarkedFrameOccu(y=y_oc, siteCovs=sc_oc)
    # summary(umf_oc0)

    # Add availability
     # Simulate survey durations for OC data
    doc <- runif(nsites_OC, min(range.dur.OC), max(range.dur.OC))

    # Linear model for the singing/activity rate phi
    phi <- mean.phi              # possibly later also add covariate

    pdoc <- 1-exp(-1*doc*phi) # Compute availability prob.
    umf_oc1 <- umf_oc0             # Copy unmarked data frame
    n_ind_avail <- rbinom(nsites_OC, dat3.raw$counts, pdoc)
    # y_oc2 <- matrix(as.numeric(n_ind_avail > 0), ncol=1
    umf_oc1@y <- matrix(as.numeric(n_ind_avail > 0), ncol=1)
  }

  # Make plot of availability model
  if(show.plots) {
    dur <- 0:max(range.dur.PC, range.dur.OC)
    par(mfrow = c(1,1), mar = c(5,5,4,3), cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
    plot(dur, 1-exp(-1*dur*mean.phi), xlab = 'Survey duration (min)', ylab = 'Availability probability',
        main = paste('Parametric function for availability,\nwith phi =', round(mean.phi,2)),
        type = 'l', frame = FALSE, ylim = c(0, 1), lwd = 5, col = 'grey')
  }

  # Numerical output
  out <- list( # input arguments
      mean.density = mean.density, beta.density = beta.density, mean.phi = mean.phi, beta.phi = beta.phi,
      nsites_HDS = nsites_HDS, sigHDS = sigHDS, nbins = nbins, maxDist_HDS = maxDist_HDS, range.dur.HDS = range.dur.HDS,
      nsites_PC = nsites_PC, sigPC = sigPC, maxDist_PC = maxDist_PC, range.dur.PC = range.dur.PC,
      nsites_OC = nsites_OC, sigOC = sigOC, maxDist_OC = maxDist_OC, range.dur.OC = range.dur.OC,
      # generated values
      phi = phi, dds = dds, umf_hds0 = umf_hds0, umf_hds1 = umf_hds1)
  if(nsites_PC > 0) {
      out$dpc <- dpc
      out$umf_pc0 <- umf_pc0
      out$umf_pc1 <- umf_pc1
  }
  if(nsites_OC > 0) {
      out$doc <- doc
      out$umf_oc0 <- umf_oc0
      out$umf_oc1 <- umf_oc1
  }

  return(out)

}  # ------------ End of function definition ---------------
