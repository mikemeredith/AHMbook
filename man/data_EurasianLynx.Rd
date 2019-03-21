\name{EurasianLynx}
\alias{EurasianLynx}
\encoding{UTF-8}
\docType{data}
\title{
Data for Eurasian lynx from Switzerland and Italy.
}
\description{
The data are observations of Eurasian lynx (\emph{Lynx lynx}) by experts and by the general public, summarised to 10 x 10 km grid and 3 occasions (Nov-Dec, Jan-Feb, Mar-Apr). The observations are classified as 'certain' or 'uncertain', the latter with false positives. This data set  covers Switzerland and Italy from 1994 to 2016.
}

\usage{data("EurasianLynx")}

\format{
\code{EurasianLynx} is a data frame with 43,332 rows corresponding to observations, and 7 columns:

\describe{
  \item{type }{factor: 'certain' or 'uncertain'.}
  \item{site.nr }{site (10 x 10 km cell) identifier.}
  \item{y.1, y.2, y.3 }{detection (1)/non-detection (0) for the three occasions.}
  \item{Year }{the year of the observation.}
  \item{Cntry }{the country, 'Italy' or 'Switzerland'.}
}
}
\source{
The data were collected by the Swiss Foundation KORA Carnivore Ecology and Wildlife Management and Progetto Lince Italia and obtained from Fridolin Zimmermann and Anja Molinari-Jobin.
}
\references{
Molinari-Jobin, A., Kéry, M., Marboutin, E., Molinari, P., Koren, I., Fuxjäger, C., Breitenmoser-Würsten, C., Wölfl, S., Fasel, M., Kos, I., Wölfl, M., & Breitenmoser, U. (2012) Monitoring in the presence of species misidentification: the case of the Eurasian lynx in the Alps. \emph{Animal Conservation}, 15, 266-273.

Kéry & Royle (2019) \emph{Applied Hierarchical Modeling in Ecology} Section 18.3.2
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
