\name{SwissTits}
\alias{SwissTits}
\encoding{UTF-8}
\docType{data}
\title{
Data from the Swiss Breeding Bird Survey for 6 species of tits
}
\description{
The Swiss breeding bird survey ("Monitoring Häufige Brutvögel" MHB) has monitored the populations of approximately 150 common species since 1999. The MHB sample consists of 267 1-km squares that are laid out as a grid across Switzerland. Fieldwork is conducted by about 200 skilled birdwatchers annually, most of them volunteers. Avian populations are monitored using a simplified territory mapping protocol (also called spot-mapping West of the Atlantic), where each square is surveyed up to three times during the breeding season (only twice above the tree line around 2000 m elevation). Surveys are conducted along a transect that does not change over the years and is 4-6 kms in length on average.

The list \code{SwissTits} has the data for six species of tits from 2004 to 2013. There are some missing values: see Details.
}

\usage{data("SwissTits")}

\format{
\code{SwissTits} is a list with 5 elements:

\describe{
  \item{\bold{species }}{a data frame with a row for each species and the following columns:

    \itemize{
      \item specid : a numeric species ID based on phylogeny.
      \item latname : a 6-letter abbreviation of the Latin name.
      \item name : the English name.
    }}}
\describe{
  \item{\bold{sites }}{a data frame with a row for each 1x1 km quadrat and the following columns:

    \itemize{
      \item siteID : an alphanumeric site identifier.
      \item coordx : the x coordinate of the center of the quadrat; the coordinate reference system intentionally not specified.
      \item coordy : the y coordinate of the center of the quadrat.
      \item AQ : an identifier for the 10km x 10km block called "Atlas Quadrat" within which the site falls.
      \item AQ.coordx : the x coordinate of the center of the AQ.
      \item AQ.coordy : the y coordinate of the center of the AQ.
      \item elev : the mean elevation of the quadrat, m.
      \item rlength : the length of the route walked in the quadrat, km.
      \item forest : percentage forest cover.
    }}}
\describe{
  \item{\bold{counts }}{a sites x replicates x years x species array of counts}

  \item{\bold{date }}{a sites x replicates x years array with Julian dates of the surveys, 1 April = 1}

  \item{\bold{dur }}{a sites x replicates x years array with the duration of each survey, mins}
}
}
\details{

Missing values in the \code{date} array indicate that the corresponding survey was not carried out.

On 26 occasions when surveys were carried out, the duration was not recorded, resulting in additional NAs in the \code{dur} array.

A new method for recording breeding territories was introduced in 2004, but the old protocol was in use at some sites until 2013. Surveys with the old protocol have the counts shown as NA in the \code{count} array.
}

\note{
Sections 6.9.1 and 6.13.1 of the AHM1 book have code to read in data from a CSV file, "SwissTits_mhb_2004_2013.csv". The \code{SwissTits} list has all the same data in a more compact format. See Examples for ways to generate the objects used in the book from the list.
}
\source{
Swiss Ornithological Institute
}
\references{
Kéry, M. & Royle, J.A. (2016) \emph{Applied Hierarchical Modeling in Ecology} AHM1 - 6.9.1 and 6.13.1.
}
\examples{
data(SwissTits)
str(SwissTits)

# Section 6.9.1, foot of p.257 and top of p.258
# -------------
y0 <- SwissTits$counts[, , '2013', 'Great tit']
( NA.sites <- which(rowSums(is.na(y0)) == 3) ) # Unsurveyed sites
y <- y0[-NA.sites, ]                 # Drop them from the count data
tits <- SwissTits$sites[-NA.sites, ] # Also drop from the site covariates
str(y)
# Get date and duration data for 2013, without the NA.sites rows:
date <- SwissTits$date[-NA.sites, , '2013']
dur <- SwissTits$dur[-NA.sites, , '2013']

# Section 6.13.1, p.303
# --------------
# Get the count data for 2013 (all species)
y0 <- SwissTits$count[, , '2013', ]
str(y0)
# We keep the sites with count data, remove those with 3 NAs
# See which sites have counts in 2013 for (say) Great tits:
keep <- which(rowSums(is.na(y0[, , "Great tit"])) != 3)
length(keep)
y <- y0[keep, , ]
# Get the covariate data for the 'keep' sites
elev <- SwissTits$sites$ele[keep]
route <- SwissTits$sites$rlength[keep]
forest <- SwissTits$sites$forest[keep]
date <- SwissTits$date[keep, , '2013']  # Survey date
dur <- SwissTits$dur[keep, , '2013']    # Survey duration

# Degrade counts to detection/nondetection data
y3DRN <- y
y3DRN[y3DRN > 1] <- 1
str(y3DRN)
# Final detail...
( spec.names <- paste0(SwissTits$species$name, "s") )
}
\keyword{datasets}
