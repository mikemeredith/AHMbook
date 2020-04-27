
# Some crude checks for the simulation functions in AHMbook package

# This uses the new improved version of checkTotal
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

# halt <- FALSE
halt <- TRUE
  # If TRUE, the script will stop if output is not correct.

data("BerneseOberland")
  # checkTotal(BerneseOberland)
  checkTotal(BerneseOberland, 236659.42128, 5, halt)

data("crestedTit")
  # checkTotal(crestedTit)
  checkTotal(crestedTit, 171.41585, 5, halt)

data("crossbillAHM")
  # checkTotal(crossbillAHM)
  checkTotal(crossbillAHM, 15801.46998, 5, halt)

data("cswa")
  # checkTotal(cswa)
  checkTotal(cswa, 969.54437, 5, halt)

#data("dragonflies")  # Doesn't work

data("duskySalamanders")
  # checkTotal(duskySalamanders)
  checkTotal(duskySalamanders, 1.32143, 5, halt)

data("EurasianLynx")
  # checkTotal(EurasianLynx)
  checkTotal(EurasianLynx, 935.49487, 5, halt)

data("Finnmark")
  # checkTotal(Finnmark)
  checkTotal(Finnmark, 163.55519, 5, halt)

data("FrenchPeregrines")
  # checkTotal(FrenchPeregrines)
  checkTotal(FrenchPeregrines, 3.11435, 5, halt)

data("greenWoodpecker")
  # checkTotal(greenWoodpecker)
  checkTotal(greenWoodpecker, 8653.59173, 5, halt)

data("HubbardBrook")
  # checkTotal(HubbardBrook)
  checkTotal(HubbardBrook, 77.45959, 5, halt)

data("MesoCarnivores")
  # checkTotal(MesoCarnivores)
  checkTotal(MesoCarnivores, 0.41693, 5, halt)

data("MHB2014")
  # checkTotal(MHB2014)
  checkTotal(MHB2014, 34974.37985, 5, halt)

data("spottedWoodpecker")
  # checkTotal(spottedWoodpecker)
  checkTotal(spottedWoodpecker, 171055.65596, 5, halt)

data("SwissAtlasHa")
  # checkTotal(SwissAtlasHa)
  checkTotal(SwissAtlasHa, 57931.50771, 5, halt)

data("SwissEagleOwls")
  # checkTotal(SwissEagleOwls)
  checkTotal(SwissEagleOwls, 706.56796, 5, halt)

data("SwissMarbledWhite")
  # checkTotal(SwissMarbledWhite)
  checkTotal(SwissMarbledWhite, 137.93353, 5, halt)

## SwissSquirrels is in a txt file
fn <- file.path(find.package("AHMbook"), "extdata", "SwissSquirrels.txt")
ss <- read.table(fn, header = TRUE)  # checkTotal(SwissTits)
  # checkTotal(ss)
  checkTotal(ss, 86839.57389, 5, halt)

data("SwissTits")
  # checkTotal(SwissTits)
  checkTotal(SwissTits, 64312.21470, 5, halt)

data("treeSparrow")
  # checkTotal(treeSparrow)
  checkTotal(treeSparrow, 88318.33976, 5, halt)

data("ttdPeregrine")
  # checkTotal(ttdPeregrine)
  checkTotal(ttdPeregrine, 48.49225, 5, halt)

data("UKmarbledWhite")
  # checkTotal(UKmarbledWhite)
  checkTotal(UKmarbledWhite, 28116.73671, 5, halt)

data("wagtail")
  # checkTotal(wagtail)
  checkTotal(wagtail, 45.89299, 5, halt)

data("wigglyLine")
  # checkTotal(wigglyLine)
  checkTotal(wigglyLine, 2.0322, 5, halt)

data("waterVoles")
  # checkTotal(waterVoles)
  checkTotal(waterVoles, 11.96832, 5, halt)

data("willowWarbler")
  # checkTotal(willowWarbler)
  checkTotal(willowWarbler, 99346.51996, 5, halt)

if(halt) cat("\n\nYay! No problems!\n\n")
