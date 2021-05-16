\name{modSelFP}
\alias{modSelFP}
\encoding{UTF-8}
\title{
Model selection results from a list of 'unmarkedFitOccuFP' models
}
\description{
The \pkg{unmarked} functions \code{\link[unmarked]{fitList}} and \code{\link[unmarked]{modSel}} do now work for 'unmarkedFitOccuFP' models, so this stop-gap function is deprecated. Use \code{modSel(fitList(fits=mod.list))} instead.
}
\usage{
modSelFP(mod.list)
}

\arguments{
  \item{mod.list}{
A list of fitted model objects of class 'unmarkedFitOccuFP' as returned by \pkg{unmarked}.
}
}

\value{
A data frame with a row for each model in the input list and columns for number of parameters, AIC, deltaAIC, AIC weight and cumulative AIC weight.
}
\references{
Kéry, M. & Royle, J.A. (2021) \emph{Applied Hierarchical Modeling in Ecology} AHM2 - 7.2.2.
}
