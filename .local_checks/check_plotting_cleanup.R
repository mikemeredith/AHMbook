
# Check that plotting functions plot and clean up pars and devAskNew

# Source the script and see if the before and after plots are identical. If not,
#   look in "Cleanup_check.pdf" for the offending function.

# REMOVE the "Cleanup_check.pdf" file before committing to Github!

library(AHMbook)

graphics.off()  # Start with a fresh window

testPlot <- function(label) {
  stdNames <- c("black", "red", "green3", "blue", "cyan",
    "magenta", "yellow",  "gray")
  plot(1:8, rep(1,8), type='n', xlim=c(0.5, 8.5),
    ylim=c(0.5, 2.5), xlab="Check margins", ylab="Check orientation",
    main="Check size and font")
  if(!missing(label))
    mtext(label)
  # colours
  text(c(4,4), c(1.1, 2.1),
    c("Colours by number (depends on palette)", "Colours by name"))
  points(1:8, rep(1,8), pch=16, col=1:8)
  segments(1:8, rep(0.7, 8), 1:8, rep(0.9, 8), col=1:8)
  points(1:8, rep(2,8), pch=16, col=stdNames[1:8])
  segments(1:8, rep(1.7, 8), 1:8, rep(1.9, 8), col=1:8)
  text(c(1.2), c(1.3, 1.5), c("color #2", "color 'red'"), pos=4)
  points(c(1,1), c(1.3, 1.5), pch=15, cex=5, col=c(2, 'red'))
  # Point and  line size/width/type
  text(4, c(1.3,1.5), c("defaults", "current settings"), pos=4)
  points(3, 1.3, pch=1, cex=1, col='black', lwd=1, lty=1)
  points(3, 1.5)
  segments(3.2, 1.5, 3.9, 1.5)
  segments(3.2, 1.3, 3.9, 1.3, cex=1, col='black', lwd=1, lty=1)

  text(4, 2.4, "Check position: __centered on dot.")
  points(4, 2.4)
  text(4, 2.3, "Check point size, line width, line type.")
  mtext("Normal margin v", side =1, line=4, adj=0)
  mtext("Normal margin ^", side =2, line=3, adj=0)
  mtext("Normal margin ^", side =3, line=3, adj=0)
  mtext("Normal margin v", side =4, line=1, adj=0)
}

testPlot("before")

pdf("Cleanup_check.pdf")
# AHM1 Chapter 1
# --------------
tmp <- sim.fn()  # 1 page
testPlot() # ok
tmp <- sim.fn(show.plot=FALSE) # ok

# AHM1 Chapter 4
# ---------
tmp <- data.fn() # 3 pages
testPlot() # ok
tmp <- data.fn(show.plot=FALSE) # ok

# AHM1 Chapter 6
# -------
tmp <- simNmix()  # 9 pages
testPlot() # ok
tmp <- simNmix(show.plot=FALSE) # ok
## warnings
tmp <- simNmix(Neg.Bin=TRUE)
tmp <- simNmix(open.N = TRUE)
tmp <- simNmix(Neg.Bin=TRUE, open.N=TRUE)
testPlot() # ok
areas <- runif(267, 1, 2) # Use same areas for all
tmp <- simNmix(area=areas)
tmp <- simNmix(area=areas, Neg.Bin=TRUE)
tmp <- simNmix(area=areas, open.N = TRUE)
tmp <- simNmix(area=areas, Neg.Bin=TRUE, open.N=TRUE)
testPlot() # ok

tmp <- simpleNmix()  # 1 page
testPlot() # ok
tmp <- simpleNmix(show.plot=FALSE)  # ok

tmp <- playRN()  # 1 page, 1 plot
testPlot() ## OK!
tmp <- playRN(show.plot=FALSE)
## warnings
# Try this too:
par(mfrow = c(2,2))
for(i in 1:4)
  tmp <- playRN()
par(mfrow = c(1,1))

# AHM1 Chapter 8
# ---------
tmp <- sim.ldata()  # 1 page
testPlot() # ok
tmp <- sim.ldata(show.plot=FALSE)  # ok

tmp <- sim.pdata()  # 1 page
testPlot() # ok
tmp <- sim.pdata(show.plot=FALSE)  # ok

tmp <- simHDS("point")  # 1 page
tmp <- simHDS("line")  # 1 page
testPlot() # ok
tmp <- simHDS(show.plot=FALSE)  # ok

# AHM1 Chapter 9
# ---------
tmp <- simHDSg("line")  # 1 page
testPlot() # ok
tmp <- simHDSg("point")  # 1 page
testPlot() # ok
tmp <- simHDSg(show.plot=FALSE)  # ok

tmp <- simHDStr("line", "removal")  # 1 page
tmp <- simHDStr("line", "double")  # 1 page
tmp <- simHDStr("point", "removal")  # 1 page
tmp <- simHDStr("point", "double")  # 1 page
testPlot() # ok
tmp <- simHDStr(show.plot=FALSE)  # ok

## tmp <- simHDSopen()  # no plots

tmp <- sim.spatialDS(model="half")  # 1 page
tmp <- sim.spatialDS(B=10, model="half")  # 1 page
tmp <- sim.spatialDS(model="logit")  # 1 page
testPlot() ## ok
tmp <- sim.spatialDS(show.plot=FALSE)  # ok

tmp <- sim.spatialHDS() # default 3 plots
tmp <- sim.spatialHDS(B=50) # default 3 plots
testPlot() ## ok
tmp <- sim.spatialHDS(show.plot=FALSE)  # ok

# AHM1 Chapter 10
# ----------
tmp <- simOcc()  # 2 pages
testPlot() # ok
tmp <- simOcc(show.plot=FALSE)  # ok

tmp <- sim3Occ()  # 1 page
testPlot() # ok
tmp <- sim3Occ(show.plot=FALSE)  # ok

tmp <- simOccttd()  # 1 page
testPlot() ## ok
tmp <- simOccttd(show.plot=FALSE)  # ok

tmp <- wigglyOcc()  # 1 page
testPlot() # ok
tmp <- wigglyOcc(show.plot=FALSE)  # ok

# AHM1 Chapter 11
# ----------
tmp <- simComm(type="counts")  # 2 pages
testPlot() # ok
tmp <- simComm(type="counts", show.plot=FALSE)  # ok
tmp <- simComm(type="det/nondet")  # 2 pages
testPlot() # ok
tmp <- simComm(type="det/nondet", show.plot=FALSE)  # ok

# AHM2 Chapter 1
# ----------
tmp <- simNpC()  # 1 page
testPlot() # ok
tmp <- simNpC(show.plot=FALSE)  # ok

tmp <- simPOP()  # 3 pages
testPlot() # ok
tmp <- simPOP(show.plot=FALSE)  # ok

tmp <- simPH()  # 4 pages
testPlot() # ok
tmp <- simPH(show.plot=FALSE)  # ok

# AHM2 AHM1 Chapter 2
# ----------
# tmp <- simDM0()  # 0 pages

tmp <- simDM()  # 1 page
testPlot() # ok
tmp <- simDM(show.plot=FALSE)  # ok

tmp <- simMultMix()  # 1 page
testPlot() # ok

# AHM2 Chapter 3
# ----------
tmp <- simCJS()  # 2 pages
testPlot() # ok
tmp <- simCJS(show.plot=FALSE)  # ok

# AHM2 Chapter 4
# ----------
tmp <- simDynocc()  # 2 pages
testPlot() # ok
tmp <- simDynocc(show.plot=FALSE)  # ok

tmp <- simDemoDynocc()  # 1 page
testPlot() # ok
tmp <- simDemoDynocc(show.plot=FALSE)  # ok

# AHM2 Chapter 5
# ----------
tmp <- simDCM()  # 2 pages
testPlot() # ok
tmp <- simDCM(show.plot=FALSE)  # ok

# AHM2 Chapter 9
# ----------
tmp <- simDynoccSpatial(nyears=3)  # 7 pages
testPlot() # ok

tmp <- simOccSpatial()  # 4 pages
testPlot() # now ok

tmp <- simNmixSpatial()  # 4 pages
testPlot() # now ok

# AHM2 Chapter 10
# ----------
tmp <- simPPe()  # 1 pages
testPlot() # now ok

tmp <- simDataDK()  # 2 pages
testPlot() # ok

# AHM2 Chapter 11
# ----------
tmp <- simSpatialDSline()  # 1 pages
testPlot() # now ok

# tmp <- sim.spatialDSte()  # no plots
# testPlot() # now ok

X <- matrix(c(0.46, 1.09, 1.72, 1.93, 2.52, 2.85, 2.19, 1.45, 1.70, 2.24,
  3.48, 3.74, 3.34, 2.66, 2.37, 1.81, 1.63, 1.54, 1.24, 0.79), ncol=2)
tmp <- simDSM(X)  # 1 pages
testPlot() # now ok

dev.off()

testPlot("after")

