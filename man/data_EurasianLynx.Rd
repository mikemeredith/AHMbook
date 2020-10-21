\name{EurasianLynx}
\alias{EurasianLynx}
\encoding{UTF-8}
\docType{data}
\title{
Data for Eurasian lynx from Switzerland and Italy (only Alps).
}
\description{
The data are observations of Eurasian lynx (\emph{Lynx lynx}) confirmed by experts and by the general public (and not confirmed by experts), summarized to a 10 x 10 km grid and for 3 occasions during a 'winter' (Nov-Dec, Jan-Feb, Mar-Apr). The observations are classified as 'certain' or 'uncertain', and only the latter are assumed to be contaminated with false positives (but both are subject to false-negative errors). This data set covers the Alps in Switzerland and Italy from 1994 to 2016.
}

\usage{data("EurasianLynx")}

\format{
\code{EurasianLynx} is a data frame with 43,332 rows corresponding to observations, and 10 columns:

\describe{
  \item{type }{factor: 'certain' or 'uncertain'.}
  \item{site.nr }{site (10 x 10 km cell) identifier.}
  \item{y.1, y.2, y.3 }{detection (1)/non-detection (0) for the three occasions.}
  \item{Year }{the year of the observation.}
  \item{Cntry }{the country, 'Italy' or 'Switzerland'.}
  \item{xcoord }{the x coordinate of the cell, km E of the origin.}
  \item{ycoord }{the y coordinate of the cell, km N of the origin.}
  \item{forest }{the percentage forest cover in the cell.}
}
}
\source{
The Foundation KORA (SCALP Project) and the Progetto Lince Italia (Anja Molinari-Jobin, Urs Breitenmoser, Fridolin Zimmermann).
}
\references{
Molinari-Jobin, A., Kéry, M., Marboutin, E., Molinari, P., Koren, I., Fuxjäger, C., Breitenmoser-Würsten, C., Wölfl, S., Fasel, M., Kos, I., Wölfl, M., & Breitenmoser, U. (2012) Monitoring in the presence of species misidentification: the case of the Eurasian lynx in the Alps. \emph{Animal Conservation}, 15, 266-273.

Kéry, M. & Royle, J.A. (2021) \emph{Applied Hierarchical Modeling in Ecology} AHM2 - 7.3.2 and 9.6.1.3.
}
\examples{
data(EurasianLynx)
lynx <- EurasianLynx # Shorter name
str(lynx)

# Create additional columns needed for the analysis in 'unmarked':
lynx$occ.1 <- 1
lynx$occ.2 <- 2
lynx$occ.3 <- 3
lynx$sYear <- standardize(lynx$Year)
str(lynx)
}
\keyword{datasets}
