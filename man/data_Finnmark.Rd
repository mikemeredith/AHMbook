\name{Finnmark}
\alias{Finnmark}
\encoding{UTF-8}
\docType{data}
\title{
Data from surveys of birds in Finnmark in NE Norway
}
\description{
A total of 37 plots were placed spanning the existing variation in extent and fragmentation of willow thickets. Each plot was visited 3-5 times in early July each year from 2005 to 2008. Birds observed within 100m during a 15 min period were recorded.
}

\usage{data("Finnmark")}

\format{
\code{Finnmark} is a list with 4 elements:

\describe{
  \item{\bold{species }}{a data frame with a row for each species and the following columns:

    \itemize{
      \item species : the English name.
      \item latin : the Latin name.
      \item assemblage : the guild to which each species belongs.
    }}}
\describe{
  \item{\bold{sites }}{a data frame with a row for each plot and the following columns:

    \itemize{
      \item region : a factor, 3 levels, Ilford, Komag, and Vestre Jakobselv.
      \item catchment : a factor, 11 levels.
      \item plot : a factor, the alphanumeric ID of each plot.
      \item plotnr : the ID number of each plot.
      \item area : the percentage of a 400 x 400m quadrat centered on the sampling point covered by tall willow thickets.
      \item edge : the length (m) of the edge of the willow thickets in the quadrat.
      \item height : the mean height of willows at 4 measuring points.
      \item density : a measure of thicket density at 4 measuring points.
    }}}
\describe{
  \item{\bold{counts }}{a sites x replicates x years x species array of counts}

  \item{\bold{timeOfDay }}{a sites x replicates x years array with the time of day that the survey was carried out.}
}
}

\source{
Data generously provided by John-André Henden.
}
\references{
Henden, J.-A., Yoccoz, N.G., Ims, R.A., Langeland, K. (2013) How Spatial Variation in Areal Extent and Configuration of Labile Vegetation States Affect the Riparian Bird Community in Arctic Tundra. PLoS ONE 8(5): e63312.

Kéry, M. & Royle, J.A. (2021) \emph{Applied Hierarchical Modeling in Ecology} AHM2 - 5.
}
%% will be able to tidy up the examples with page numbers etc once we have the proofs.
\examples{
data(Finnmark)
str(Finnmark)

# Create the objects needed for the analysis:
# Extract the numeric covariates and standardise
str(Finnmark$sites)
scov <- scale(Finnmark$sites[, 5:8])
str(scov)

# Convert the 'counts' array to detection/nondetection data:
y <- Finnmark$counts > 0
storage.mode(y) <- "integer"
str(y)

# Get the guild information for each species
guild <- Finnmark$species$assemblage

# Standardise the timeOfDay and replace NAs with 0
tod <- with(Finnmark, (timeOfDay - mean(timeOfDay, na.rm=TRUE))/ sd(timeOfDay, na.rm=TRUE))
tod[is.na(tod)] <- 0
str(tod)

}
\keyword{datasets}
