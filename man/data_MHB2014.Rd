\name{MHB2014}
\alias{MHB2014}
\encoding{UTF-8}
\docType{data}
\title{
Data from the Swiss Breeding Bird Survey 2014
}
\description{
The Swiss breeding bird survey ("Monitoring Häufige Brutvögel" MHB) has monitored the populations of 150 common species since 1999. The MHB sample consists of 267 1-km squares that are laid out as a grid across Switzerland. Fieldwork is conducted by about 200 skilled birdwatchers, most of them volunteers. Avian populations are monitored using a simplified territory mapping protocol, where each square is surveyed up to three times during the breeding season (only twice above the tree line). Surveys are conducted along a transect that does not change over the years.

The list \code{MHB2014} has the data for 2014.
}

\usage{data("MHB2014")}

\format{
\code{MHB2014} is a list with 5 elements:

\describe{
  \item{\bold{species }}{a data frame with row for 158 species, including 13 species not recorded in 2014, and the following columns:

    \itemize{
      \item specid : a numeric species ID based on phylogeny.
      \item latabb : a 6-letter abbreviation of the Latin name.
      \item engname : the English name.
      \item latname : the Latin name.
      \item body.length : body length in cm.
      \item body.mass : body mass in g.
      \item wing.span : wing span in cm.
    }}}
\describe{
  \item{\bold{sites }}{a data frame with rows for 267 1x1 km quadrat, including 1 quadrat not surveyed in 2014, and the following columns:

    \itemize{
      \item siteID : an alphanumeric site identifier.
      \item coordx : the x coordinate of the centre of the quadrat; the coordinate reference system intentionally not specified.
      \item coordy : the y coordinate of the centre of the quadrat.
      \item elev : the mean elevation of the quadrat, m.
      \item rlength : the length of the route walked in the quadrat, km.
      \item nsurvey : the number of replicate surveys \emph{planned} in the quadrat; above the tree-line 2, otherwise 3.
      \item forest : percentage forest cover.
      \item obs14 : identifying number of the observer.
    }}}
\describe{
  \item{\bold{counts }}{a sites x replicates x species array of counts}

  \item{\bold{date }}{a sites x replicates matrix with Julian dates of the surveys, 1 April = 1}

  \item{\bold{dur }}{a sites x replicates matrix with the duration of each survey, mins}
}
}

\note{
Section 11.3 of the book has code to read in data from a CSV file, "MHB_2014.csv". This is a huge file, because the site data are repeated for all 158 species and the species data are repeated for all 267 sites. The \code{MHB2014} list has all the same data in a more compact format. See Examples for ways to generate the objects used in the book from the list.
}
\source{
Swiss Ornithological Institute
}
\references{
Kéry & Royle (2016) \emph{Applied Hierarchical Modeling in Ecology} Section 11.3
}
\examples{
data(MHB2014)
str(MHB2014)

# Create the objects at foot of p.644:
( nsite <- nrow(MHB2014$sites) )    # number of sites in Swiss MHB
nrep <- 3                           # maximum number of replicate surveys per season
( nspec <- nrow(MHB2014$species) )  # 158 species occur in the 2014 data
# Check the dimensions of the 'count' array:
dim(MHB2014$count) == c(nsite, nrep, nspec)

# Create the detection/nondetection matrix 'y':
y <- MHB2014$count > 0        # this is logical, convert to integer
storage.mode(y) <- "integer"  # don't use 'as.integer', that strips out dimensions and names
str(y)

# Pull out and check the data for common chaffinch, p.645:
head(tmp <- y[, , "Common Chaffinch"])
tail(tmp)

}
\keyword{datasets}
