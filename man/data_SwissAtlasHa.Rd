\name{SwissAtlasHa}
\alias{SwissAtlasHa}
\encoding{UTF-8}
\docType{data}
\title{
Data for 1ha quadrats from the Swiss Breeding Bird Atlas
}
\description{
A 1ha-scale subset of the count data from the most recent Swiss Breeding Bird Atlas (Knaus et al. 2018). The data consists of replicate counts of birds from 2,318 1ha quadrats for 78 passerine species. Each quadrat is surveyed up to three times during the breeding season (only twice above the tree line).
}

\usage{data("SwissAtlasHa")}

\format{
\code{SwissAtlasHa} is a list with 3 elements:

\describe{
  \item{\bold{counts }}{a sites x replicates x species array of counts}

  \item{\bold{sites }}{a data frame with rows for 2318 1ha quadrats and the following columns:

    \itemize{
       \item x : the x coordinate of the center of the quadrat; the coordinate reference system is intentionally not specified.
      \item y : the y coordinate of the center of the quadrat.
      \item elevation : the mean elevation of the quadrat, m.
      \item slope : the slope.
      \item north : the aspect, with due south = -1, due north = 1.
      \item forest : proportion of forest cover.
      \item nsurvey : the number of replicate surveys carried out in the quadrat.
    }}

  \item{\bold{date }}{a sites x replicates matrix with Julian dates of the surveys, 1 January = 1}
}
}

\source{
Swiss Ornithological Institute
}
\references{
Knaus, P., S. Antoniazza, S. Wechsler, J. Guélat, M. Kéry, N. Strebel, & T. Sattler (2018) Brutvogelatlas 2013–2016. Bestandsentwicklung der Brutvögel der Schweiz und des Fürstentums Liechtensteins (Swiss Breeding Bird Atlas 2013–2016). Schweizerische Vogelwarte, Sempach.

Kéry, M. & Royle, J.A. (2021) \emph{Applied Hierarchical Modeling in Ecology} AHM2 - 8.4.2.
}
\examples{
data(SwissAtlasHa)
str(SwissAtlasHa)
# display the species names:
dimnames(SwissAtlasHa$counts)[[3]]
}
\keyword{datasets}
