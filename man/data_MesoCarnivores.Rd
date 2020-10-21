\name{MesoCarnivores}
\alias{MesoCarnivores}
\encoding{UTF-8}
\docType{data}
\title{
Camera trap data for Bobcat, Coyote and Red Fox
}
\description{
Camera trap (detection/nondetection) data for Bobcat (\emph{Lynx rufus}), Coyote (\emph{Canis latrans}) and Red Fox (\emph{Vulpes vulpes}) from six Mid-Atlantic states in the eastern United States analyzed by Rota et al (2016).
}

\usage{data("MesoCarnivores")}

\format{
\code{MesoCarnivores} is a list with 4 elements:

\describe{
  \item{\bold{bobcat, coyote, redfox }}{1/0 detection data for the respective species: matrices with rows for 1437 sites x 3 replicates, where each replicate corresponds to 1 week of observations.}

  \item{\bold{sitecov }}{a data frame with rows for 1437 sites and the following columns:

    \itemize{
      \item Dist_5km : the proportion of disturbed land in the surrounding 5km.
      \item HDens_5km : housing density in the surrounding 5km.
      \item Latitude : the latitude of the point, decimal degrees divided by 100.
      \item Longitude : the longitude, decimal degrees divided by 100.
      \item People_site : number of people photographed at a site divided by 1000.
      \item Trail : 1 if the camera was located on a trail, 0 otherwise.
    }}
}
}

\source{
Data courtesy of eMammal, Roland Kays, Arielle Parsons, and their group at the North Carolina Museum of Natural Sciences and North Carolina State University.
}
\references{
Rota, C.T., Ferreira, M.A.R., Kays, R.W., Forrester, T.D., Kalies, E.L., McShea, W.J., Parsons, A.W., & Millspaugh, J.J. (2016) A multispecies occupancy model for two or more interacting species. \emph{Methods in Ecology and Evolution}, 7, 1164-1173.

Kéry, M. & Royle, J.A. (2021) \emph{Applied Hierarchical Modeling in Ecology} AHM2 - 8.2.
}
\examples{
data(MesoCarnivores)
str(MesoCarnivores)
}
\keyword{datasets}
