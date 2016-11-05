\name{jay}
\alias{jay}
\docType{data}
\title{
European jay data from the Swiss Breeding Bird Survey
}
\description{
The \pkg{unmarked} package has a very similar data set \emph{with the same name}, and it lazy-loads without needing to use \code{date(jay)}.
}
\usage{data(jay)}
\format{
  A list with 5 elements: 

  \bold{$covinfo :'data.frame':	238 obs. of  10 variables}
  
  Covariate values for the quadrats surveyed:
  
  \describe{
    \item{\code{x}}{the x coordinate (Easting) of the quadrat}
    \item{\code{y}}{the y coordinate (Northing) of the quadrat}
    \item{\code{elev}}{average elevation of the quadrat in metres}
    \item{\code{length}}{length of the sampled route through the quadrat in kilometres}
    \item{\code{durX}}{duration of the survey in minutes}
  }
  
  \bold{$capthist :'data.frame':	238 obs. of  10 variables}
  
Capture histories for quadrats with 3 surveys (\code{100, 010, 001, 110, 101, 011, 111}) or with 2 surveys (\code{10x, 01x, 11x}). Either the former or the latter are filled with NAs.

  \bold{$gridinfo :'data.frame':	41365 obs. of  7 variables}
  
  Not used in the book.
  
  \bold{$sitecovs :'data.frame':	238 obs. of  3 variables}
  
  Covariate values for the quadrats surveyed:
  
  \describe{
    \item{\code{elev}}{average elevation of the quadrat in metres}
    \item{\code{length}}{length of the sampled route through the quadrat in kilometres}
    \item{\code{forest}}{forest cover}
  }
}
\examples{
data(jay)
## maybe str(jay) ; plot(jay) ...
}
\keyword{datasets}
