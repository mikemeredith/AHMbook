\name{modSelFP}
\alias{modSelFP}
\encoding{UTF-8}
\title{
Model selection results from a list of 'unmarkedFitOccuFP' models
}
\description{
The \pkg{unmarked} functions \code{\link[unmarked]{fitList}} and \code{\link[unmarked]{modSel}} do not currently work for 'unmarkedFitOccuFP' models, so this is a stop-gap function for use until this capability is added. It will be removed from \pkg{AHMbook} once it is included in \pkg{unmarked}.
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
Kéry & Royle (2016) \emph{Applied Hierarchical Modeling in Ecology} Section 18.2.2
}
\author{
Chris Sutherland.
}


