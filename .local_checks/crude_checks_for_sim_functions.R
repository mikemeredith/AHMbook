
# Some crude checks for the simulation functions in AHMbook package

# No text should appear in the Console.

# This checks both old and new (> 3.6.0) sample algorithms

# It uses the new improved version of checkTotal
checkTotal <- function(x, target=NA, digits=3, halt=FALSE) {
  addItUp <- function(x) {
    if(is.list(x))
      x <- sapply(x, addItUp)
    x1 <- suppressWarnings(try(as.numeric(x), silent=TRUE))
    if(inherits(x1, "try-error"))
      return(0)
    x2 <- mean(x1[is.finite(x1)], na.rm=TRUE)
  }
  added <- addItUp(x)
  if(is.na(target))
    return(added)
  ok <- isTRUE(all.equal(round(target, digits), round(added, digits)))
  if(halt && !ok)
    stop("Check failed for '", deparse(substitute(x)), "'", call.=FALSE)
  ok
}


library(AHMbook)

nreps <- (getRversion() >= "3.6.0") + 1
sample.kind <- c("Rounding", "Rejection")

RNGkind("Mersenne-Twister", "Inversion")

# halt <- FALSE
halt <- TRUE
  # If TRUE, the script will stop if output is not correct.

for(i in 1:nreps) {
  if(getRversion() >= "3.6.0")
    RNGkind(sample.kind=sample.kind[i])
  cat("\n\n*** Checking with sample.kind =", sample.kind[i], "***\n\n") ; flush.console()

  # AHM1 Chapter 1
  # --------------
## sim.fn ##
  set.seed(123)
  # checkTotal(sim.fn(show.plot=FALSE))
  checkTotal(sim.fn(show.plot=FALSE), 25.24703, 5, halt)

  # AHM1 Chapter 4
  # --------------
## data.fn ##
  set.seed(123)
  # checkTotal(data.fn(show.plot=FALSE))
  checkTotal(data.fn(show.plot=FALSE), 96.5407, 4, halt)

  # AHM1 Chapter 6
  # --------------
## simNmix ##
  set.seed(123)
  # checkTotal(simNmix(show.plot=FALSE))
  checkTotal(simNmix(show.plot=FALSE, verbose=FALSE), 49.1620, 4, halt)

  set.seed(123)
  # checkTotal(simNmix(Neg.Bin=TRUE, show.plot=FALSE))
  checkTotal(simNmix(Neg.Bin=TRUE, show.plot=FALSE, verbose=FALSE), 48.8637, 4, halt)

  set.seed(123)
  # checkTotal(simNmix(open.N = TRUE, show.plot=FALSE))
  checkTotal(simNmix(open.N = TRUE, show.plot=FALSE, verbose=FALSE), 50.1045, 4, halt)

  set.seed(123)
  # checkTotal(simNmix(Neg.Bin=TRUE, open.N=TRUE, show.plot=FALSE))
  checkTotal(simNmix(Neg.Bin=TRUE, open.N=TRUE, show.plot=FALSE, verbose=FALSE),
    49.3230, 4, halt)

  set.seed(456) ; areas <- runif(267, 1, 2) # Use same areas for all
  set.seed(123)
  # checkTotal(simNmix(area=areas, show.plot=FALSE))
  checkTotal(simNmix(area=areas, show.plot=FALSE, verbose=FALSE), 61.0666, halt)

  set.seed(123)
  # checkTotal(simNmix(area=areas, Neg.Bin=TRUE, show.plot=FALSE))
  checkTotal(simNmix(area=areas, Neg.Bin=TRUE, show.plot=FALSE, verbose=FALSE),
    60.4715, 4, halt)

  set.seed(123)
  # checkTotal(simNmix(area=areas, open.N = TRUE, show.plot=FALSE))
  checkTotal(simNmix(area=areas, open.N = TRUE, show.plot=FALSE, verbose=FALSE),
    62.3927, 4, halt)

  set.seed(123)
  # checkTotal(simNmix(area=areas, Neg.Bin=TRUE, open.N=TRUE, show.plot=FALSE))
  checkTotal(simNmix(area=areas, Neg.Bin=TRUE, open.N=TRUE, show.plot=FALSE, verbose=FALSE),
    61.9605, 4, halt)


## simpleNmix ##
  set.seed(123)
  # checkTotal(simpleNmix(show.plot=FALSE))
  checkTotal(simpleNmix(show.plot=FALSE), 3.98759, 5, halt)

## playRN ##
  set.seed(123)
  # checkTotal(playRN(show.plot=FALSE))
  checkTotal(playRN(show.plot=FALSE, verbose=FALSE), 67.54638, 5, halt)

  # AHM1 Chapter 8
  # --------------
## sim.ldata, sim.pdata ##
  set.seed(123)
  # checkTotal(sim.ldata(show.plot=FALSE))
  checkTotal(sim.ldata(show.plot=FALSE), 63.63492, 5, halt)

  set.seed(123)
  # checkTotal(sim.pdata(show.plot=FALSE))
  checkTotal(sim.pdata(show.plot=FALSE), 226.4170, 4, halt)

  set.seed(123)
  # checkTotal(sim.pdata(keep.all=TRUE, show.plot=FALSE))
  checkTotal(sim.pdata(keep.all=TRUE, show.plot=FALSE), 226.5527, 4, halt)

## simHDS ##
  set.seed(123)
  # checkTotal(simHDS("point", show.plot=FALSE))
  checkTotal(simHDS("point", show.plot=FALSE), 11.50938, 5, halt)

  set.seed(123)
  # checkTotal(simHDS("line", show.plot=FALSE))
  checkTotal(simHDS("line", show.plot=FALSE), 12.09830, 5, halt)

  set.seed(123)
  # checkTotal(simHDS("point", discard0=FALSE, show.plot=FALSE))
  checkTotal(simHDS("point", discard0=FALSE, show.plot=FALSE),
    11.86080, 5, halt)

  set.seed(123)
  # checkTotal(simHDS("line", discard0=FALSE, show.plot=FALSE))
  checkTotal(simHDS("line", discard0=FALSE, show.plot=FALSE),
    12.26824, 5, halt)

  # AHM1 Chapter 9
  # --------------
## simHDSg ##
  set.seed(123)
  # checkTotal(simHDSg("line", show.plot=FALSE))
  checkTotal(simHDSg("line", show.plot=FALSE), 10.77344, 5, halt)

  set.seed(123)
  # checkTotal(simHDSg("point", show.plot=FALSE))
  checkTotal(simHDSg("point", show.plot=FALSE), 10.43742, 5, halt)

  set.seed(123)
  # checkTotal(simHDSg("line", discard0=FALSE, show.plot=FALSE))
  checkTotal(simHDSg("line", discard0=FALSE, show.plot=FALSE),
    10.91265, 5, halt)

  set.seed(123)
  # checkTotal(simHDSg("point", discard0=FALSE, show.plot=FALSE))
  checkTotal(simHDSg("point", discard0=FALSE, show.plot=FALSE),
    10.722004, 6, halt)

## simHDStr ##
  set.seed(123)
  # checkTotal(simHDStr("line", "removal", show.plot=FALSE))
  checkTotal(simHDStr("line", "removal", show.plot=FALSE),
    16.286182, 6, halt)

  set.seed(123)
  # checkTotal(simHDStr("point", "removal", show.plot=FALSE))
  checkTotal(simHDStr("point", "removal", show.plot=FALSE),
    16.420240, 6, halt)

  set.seed(123)
  # checkTotal(simHDStr("line", "double", show.plot=FALSE))
  checkTotal(simHDStr("line", "double", show.plot=FALSE),
    16.254108, 6, halt)

  set.seed(123)
  # checkTotal(simHDStr("point", "double", show.plot=FALSE))
  checkTotal(simHDStr("point", "double", show.plot=FALSE),
    16.316559, 6, halt)

  set.seed(123)
  # checkTotal(simHDStr("line", "removal", discard0=TRUE, show.plot=FALSE))
  checkTotal(simHDStr("line", "removal", discard0=TRUE, show.plot=FALSE),
    15.845868, 6, halt)

  set.seed(123)
  # checkTotal(simHDStr("point", "removal", discard0=TRUE, show.plot=FALSE))
  checkTotal(simHDStr("point", "removal", discard0=TRUE, show.plot=FALSE),
    15.545217, 6, halt)

  set.seed(123)
  # checkTotal(simHDStr("line", "double", discard0=TRUE, show.plot=FALSE))
  checkTotal(simHDStr("line", "double", discard0=TRUE, show.plot=FALSE),
    15.863605, 6, halt)

  set.seed(123)
  # checkTotal(simHDStr("point", "double", discard0=TRUE, show.plot=FALSE))
  checkTotal(simHDStr("point", "double", discard0=TRUE, show.plot=FALSE),
    15.446613, 6, halt)

## simHDSopen ##
  set.seed(123)
  # checkTotal(simHDSopen("line", discard0=TRUE))  # no plots
  checkTotal(simHDSopen("line", discard0=TRUE), 7.980065, 6, halt)

  set.seed(123)
  # checkTotal(simHDSopen("line", discard0=FALSE))
  checkTotal(simHDSopen("line", discard0=FALSE), 8.513981, 6, halt)

  set.seed(123)
  # checkTotal(simHDSopen("point", discard0=TRUE))
  checkTotal(simHDSopen("point", discard0=TRUE), 7.705087, 6, halt)

  set.seed(123)
  # checkTotal(simHDSopen("point", discard0=FALSE))
  checkTotal(simHDSopen("point", discard0=FALSE), 8.503422, 6, halt)


## sim.spatialDS ##
  set.seed(123)
  # checkTotal(sim.spatialDS(model="half", show.plot=FALSE)) #####
  chk <- c(196.79522, 201.004393)
  checkTotal(sim.spatialDS(model="half", show.plot=FALSE),
    chk[i], 6, halt)

  set.seed(123)
  # checkTotal(sim.spatialDS(model="logit", show.plot=FALSE))
  chk <- c(197.024388, 200.501680)
  checkTotal(sim.spatialDS(model="logit", show.plot=FALSE),
    chk[i], 6, halt)

  set.seed(123)
  # checkTotal(sim.spatialDS(keep.all=TRUE, model="half", show.plot=FALSE))
  chk <- c(199.540896, 200.774699)
  checkTotal(sim.spatialDS(keep.all=TRUE, model="half", show.plot=FALSE),
    chk[i], 6, halt)

  set.seed(123)
  # checkTotal(sim.spatialDS(keep.all=TRUE, model="logit", show.plot=FALSE))
  chk <- c(199.546714, 200.779699)
  checkTotal(sim.spatialDS(keep.all=TRUE, model="logit", show.plot=FALSE),
    chk[i], 6, halt)

## sim.spatialHDS ##
  set.seed(123)
  # checkTotal(sim.spatialHDS( show.plot=FALSE))
  chk <- c(20.790359, 20.754315)
  checkTotal(sim.spatialHDS( show.plot=FALSE),
    chk[i], 6, halt)

  # AHM1 Chapter 10
  # ----------
## simOcc ##
  set.seed(123)
  # checkTotal(simOcc(show.plot=FALSE))
  checkTotal(simOcc(show.plot=FALSE), 17.939149, 6, halt)

## sim3Occ ##
  set.seed(123)
  # checkTotal(sim3Occ(show.plot=FALSE))
  checkTotal(sim3Occ(show.plot=FALSE, verbose=FALSE), 12.235466, 6, halt)

## simOccttd ##
  set.seed(123)
  # checkTotal(simOccttd(show.plot=FALSE))
  checkTotal(simOccttd(show.plot=FALSE, verbose=FALSE), 28.646766, 6, halt)

## wigglyOcc ##
  set.seed(123)
  # checkTotal(wigglyOcc(show.plot=FALSE))
  chk <- c(60.552083, 60.554444)
  checkTotal(wigglyOcc(show.plot=FALSE, verbose=FALSE), chk[i], 6, halt)

  # AHM1 Chapter 11
  # ---------------
## simComm ##
  set.seed(123)
  # checkTotal(simComm(type="count", show.plot=FALSE))
  checkTotal(simComm(type="count", show.plot=FALSE), 18.111469, 6, halt)

  set.seed(123)
  # checkTotal(simComm(type="det", show.plot=FALSE))
  checkTotal(simComm(type="det", show.plot=FALSE), 14.955703, 6, halt)

  # AHM2 Chapter 1
  # --------------
## simNpC ##
  set.seed(123)
  # checkTotal(simNpC(show.plot=FALSE))
  checkTotal(simNpC(show.plot=FALSE), 46.778571, 6, halt)

## simPOP ##
  set.seed(123)
  # checkTotal(simPOP(show.plot=FALSE))
  checkTotal(simPOP(show.plot=FALSE), 23.495983, 6, halt)

## simPH ##
  set.seed(123)
  # checkTotal(simPH(show.plot=FALSE))
  checkTotal(simPH(show.plot=FALSE), 45.641494, 6, halt)

  # AHM2 Chapter 2
  # --------------
## simDM0, simDM ##
  set.seed(123)
  # checkTotal(simDM0())
  checkTotal(simDM0(), 7.017861, 6, halt)

  set.seed(123)
  # checkTotal(simDM(show.plot=FALSE))
  checkTotal(simDM(show.plot=FALSE), 4.545674, 6, halt)

## simMultMix ##
  set.seed(123)
  # checkTotal(simMultMix())
  checkTotal(simMultMix(), 11.5825, 6, halt)

## simFrogDisease ##
  set.seed(123)
  # checkTotal(simFrogDisease())
  checkTotal(simFrogDisease(), 7.122865, 6, halt)

  # AHM2 Chapter 3
  # --------------
## simCJS ##
  set.seed(123)
  # checkTotal(simCJS(show.plot=FALSE))
  checkTotal(simCJS(show.plot=FALSE), 19.146852, 6, halt)

  # AHM2 Chapter 4
  # --------------
## simDynocc ##
  set.seed(123)
  # checkTotal(simDynocc(show.plot=FALSE))
  checkTotal(simDynocc(show.plot=FALSE), 16.369749, 6, halt)

## simDemoDynocc ##
  set.seed(123)
  # checkTotal(simDemoDynocc(show.plot=FALSE))
  checkTotal(simDemoDynocc(show.plot=FALSE), 14.413033, 6, halt)

  # AHM2 Chapter 5
  # ----------
## sim DCM ##
  set.seed(123)
  # checkTotal(simDCM(show.plot=FALSE))
  checkTotal(simDCM(show.plot=FALSE, verbose=FALSE), 5.866008, 6, halt)

  # AHM2 Chapter 9
  # ----------
## simDynoccSpatial ##
  # set.seed(123)
  # checkTotal(simDynoccSpatial(seed.XAC = 88, seed = 123, ask.plot=FALSE))
  checkTotal(suppressMessages(
      simDynoccSpatial(seed.XAC = 88, seed = 123, ask.plot=FALSE, verbose=FALSE)),
      208.8291, 6, halt)

## simOccSpatial ##
  set.seed(123)
  # checkTotal(simOccSpatial(show.plot=FALSE))
  chk <- c(40899.993733, 40900.162859)
  checkTotal(simOccSpatial(show.plot=FALSE, verbose=FALSE), chk[i], 6, halt)

## simNmixSpatial ##
  set.seed(123)
  # checkTotal(simNmixSpatial(show.plot=FALSE))
  chk <- c(45859.301649, 45859.490752)
  checkTotal(simNmixSpatial(show.plot=FALSE, verbose=FALSE), chk[i], 6, halt)

  # AHM2 Chapter 10
  # ----------
## simPPe ##
  set.seed(123)
  # checkTotal(simPPe())
  chk <- c(511.315739, 507.411458)
  checkTotal(simPPe(), chk[i], 6, halt)

## simDataDK ##
  set.seed(123)
  # checkTotal(simDataDK(show.plot=FALSE))
  chk <- c(850.123078, 851.718848)
  checkTotal(simDataDK(show.plot=FALSE), chk[i], 6, halt)

  set.seed(123)
  # checkTotal(simDataDK1(show.plot=FALSE))
  chk <- c(850.029055, 866.670939)
  checkTotal(simDataDK1(show.plot=FALSE), chk[i], 6, halt)

  # AHM2 Chapter 11
  # ---------------
## simSpatialDSline ##
  set.seed(123)
  # checkTotal(simSpatialDSline())
  chk <- c(77.558971, 77.268455)
  checkTotal(simSpatialDSline(), chk[i], 6, halt)

## simSpatialDSte ##
  set.seed(123)
  # checkTotal(simSpatialDSte(show.plot=FALSE))
  checkTotal(simSpatialDSte(show.plot=FALSE), 10.925548, 6, halt)

## simDSM ##
  X <- matrix(c(0.46, 1.09, 1.72, 1.93, 2.52, 2.85, 2.19, 1.45, 1.70, 2.24,
    3.48, 3.74, 3.34, 2.66, 2.37, 1.81, 1.63, 1.54, 1.24, 0.79), ncol=2)
  set.seed(123)
  # checkTotal(simDSM(X, show.plot=FALSE))
  chk <- c(99.5555, 99.3373)
  checkTotal(simDSM(X, show.plot=FALSE), chk[i], 4, halt)
}

if(halt) cat("\n\n*** Yay! No problems! ***\n\n")
