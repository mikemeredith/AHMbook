
setwd("D:/Github/AHMbook_package") # my desktop
setwd("C:/Github/AHMbook_package") # my laptop
dir()

# Spelling check
library(spelling)
update_wordlist(pkg = "AHMbook", confirm = TRUE)
out <- spell_check_package(pkg = "AHMbook")

# Misc checks
library(tools)
checkTnF("AHMbook")

# Install dependencies
install.packages(c("plotrix", "raster", "RandomFields", "coda",
    "unmarked", "mvtnorm", "spdep"))

# For package development
system("R CMD INSTALL AHMbook") # Use this for a "dev" install.
devtools::load_all("AHMbook")

# Create the AHMbook package
# ==========================
unlink(list.files(pattern="Rplots.pdf", recursive=TRUE))
system("R CMD build AHMbook")  # Produces the .tar.gz
pkg <- "AHMbook_0.2.3.9000.tar.gz"  # <-- fix version number here ################

# Pick one to check:
## on desktop
system(paste("R CMD check ", pkg))
system(paste("R CMD check ", pkg, "--as-cran"))  # as-cran now runs donttest
## on laptop
system(paste("R CMD check ", pkg, "--no-manual"))
system(paste("R CMD check ", pkg, "--as-cran --no-manual"))

# Pick one to install
system(paste("R CMD INSTALL ", pkg))            # install only, use on R devel
system(paste("R CMD INSTALL ", pkg, "--build")) # install and produce the .zip binary


# Try it out:
library(AHMbook)
?AHMbook


# Check graphics:
tmp <- simNmix()
data <-  simNmix(nsites = 130, nvisits = 2, mean.lam = 120, mean.p = 0.7,
            sigma.lam = 0.6, sigma.p.site = 0.26,  sigma.p.survey = 0.37,
            show.plot = TRUE)
str(data)


example(simDSM)
example(simExpCorrRF)
example(simOccSpatial)
example(BerneseOberland)

data(MesoCarnivores)
str(MesoCarnivores)
# display the species names:
data(SwissAtlasHa)
dimnames(SwissAtlasHa$counts)[[3]]

str(simPH())

library(AHMbook)
data(wigglyLine)
points <- sp::SpatialPoints( wigglyLine )
sLine <- sp::Line(points)
regpoints <- sp::spsample(sLine, 100, type = "regular")
set.seed(2027, kind = "Mersenne-Twister") # Fig 11.15 in the book
tmp <- simDSM(X = regpoints@coords)
str(tmp)

