\name{simDynoccSpatialData}
\alias{simDynoccSpatialData}
\docType{data}
\title{
Simulated data from \code{simDynoccSpatial} for AHM2 section 9.6.1.2 (page 564)
}
\description{
The data analyzed in Kéry & Royle (2021) section 9.6.1.2 was generated using \code{simDynoccSpatial} with the \pkg{RandomFields} package. That data set cannot be generated if  \pkg{RandomFields} is not available. This provides a data set to use with the code in that section.
}
\usage{data(simDynoccSpatialData)}
\format{
  A list with 45 elements. See the Value section in \code{\link{simDynoccSpatial}} for details.
}
\details{
The data set was generated on a machine with \pkg{RandomFields} installed, using the following code:

\code{simDynoccSpatialData <- simDynoccSpatial(side = 30, nyears = 10, nsurveys = 3, mean.psi1 = 0.1,
  range.phi = c(0.5, 0.5), range.gamma = c(0.2, 0.2), range.p = c(0.4, 0.4),
  beta.Xautolog = c(1, 1), seed.XAC = 1, seed = 24, ask.plot = TRUE)}
}
\source{
Simulated data.
}
\references{
Kéry, M. & Royle, J.A. (2021) \emph{Applied Hierarchical Modeling in Ecology} AHM2 section 9.6.1.2.
}
\examples{
data(simDynoccSpatialData)
dat <- simDynoccSpatialData
str(dat)
}
\keyword{datasets}
