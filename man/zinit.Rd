\name{zinit}
\alias{zinit}
\encoding{UTF-8}

\title{
Starting values for survival analysis in JAGS/BUGS
}
\description{
Generates a matrix for use as starting values for the survival indicator in a JAGS or BUGS model for survival analysis, traditionally designated \code{z}.
}
\usage{
zinit(CH)
}
\arguments{
  \item{CH}{
An individuals x time matrix of capture records, 1 if captured, 0 otherwise, no missing values.
}
}
\value{
An individuals x time matrix with 1 in each row after the first capture; all other elements NA, including the first capture occasion.
}
\references{
Kéry, M. & Royle, J.A. (2021) \emph{Applied Hierarchical Modeling in Ecology} AHM2 - 3.2.3.
}
\author{
Marc Kéry & Andy Royle
}

\examples{
# Generate a fake capture history
( ch <- matrix(rbinom(30, 1, 0.5), 6, 5) )
zinit(ch)
}
