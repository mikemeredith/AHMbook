\name{ch2marray}
\alias{ch2marray}
\encoding{UTF-8}

\title{
Converts capture histories to an m-array
}
\description{
Converts capture histories to an m-array.
}
\usage{
ch2marray(CH)
}
\arguments{
  \item{CH}{
An individuals x time matrix of capture records, 1 if captured, 0 otherwise, no missing values.
}
}
\value{
An m-array, a (years-1) x years matrix, where element [i, j] contains the number of individuals released in year i and recaptured in year j+1 (no recaptures occur in year 1). The last column gives the number of individuals released in year i and never recaptured.
}
\references{
Kéry & Royle (2019) \emph{Applied Hierarchical Modeling in Ecology} Section 14.4.1

Kéry and Schaub (2012) \emph{Bayesian population analysis using WinBUGS - a hierarchical perspective}, Academic Press.
}
\author{
Modified from Kéry and Schaub (2012).
}

\examples{
data(willowWarbler)
ch <- as.matrix(willowWarbler$birds[ , 1:11]) # extract capture histories.
dim(ch)
ch2marray(ch)
}
