\name{piFuns}
\alias{piFuns}
\alias{instRemPiFun}
\alias{crPiFun}
\alias{crPiFun.Mb}
\alias{MhPiFun}
\encoding{UTF-8}
\title{
Compute multinomial cell probabilities
}
\description{
Compute the cell probabilities used in the multinomial-Poisson models \code{multinomPois} and \code{gmultmix} in package \pkg{unmarked}. See \link[unmarked]{piFuns}.
}
\usage{
instRemPiFun(p)

crPiFun(p)

crPiFun.Mb(p)

MhPiFun(p)
}
\arguments{
  \item{p}{
matrix of detection probabilities at each site for each observation
}
}
\details{
\code{instRemPiFun} defines the relationship between the multinomial cell probabilities and the underlying detection probability parameters (i.e., a pi function) in a removal design with unequal length of sampling periods.

\code{crPiFun} defines a pi function for capture-recapture design with 3 surveys (defined in Section 7.8.5). NOTE that this is not the same as the custom \code{crPiFun} defined in Section 7.9.1.

\code{crPiFun.Mb} defines a pi function for capture-recapture design with 3 surveys  and behavioural response (introduced in Section 7.8.2)

\code{MhPiFun} defines a pi function for model with individual detection heterogeneity (introduced in Section 7.8.3)

}
\value{
TODO
}
\references{
Kéry & Royle (2016) \emph{Applied Hierarchical Modeling in Ecology} Chapter 7
}
\author{
Kéry & Royle
}

\examples{
# TODO
}
