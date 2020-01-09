\name{bigCrossCorr}
\alias{bigCrossCorr}
\encoding{UTF-8}
\title{
Cross correlations for MCMC output
}
\description{
A wrapper for \code{coda::crosscorr}, which calculates cross-correlations between posterior draws of parameters in Markov chain Monte Carlo output. When the output has hundreds of parameters, the matrix produced by \code{crosscorr} is unwieldy, and \code{bigCrossCorr} extracts only those greater than a given threshold.
}
\usage{
bigCrossCorr(x, big = 0.6, digits = 3)
}
\arguments{
  \item{x}{
an \code{mcmc} or \code{mcmc.list} object.
}
  \item{big}{
only values below -big or above +big will be returned
}
  \item{digits}{
the number of decimal places to return
}
}
\value{
A data frame with 2 columns for the names of parameters and a 3rd column with the cross-correlation.
}
\author{
Mike Meredith
}

\seealso{
\code{\link[coda]{crosscorr}} in package \pkg{coda}.
}
