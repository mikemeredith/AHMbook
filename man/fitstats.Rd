\name{fitstats}
\alias{fitstats}
\alias{fitstats2}
\encoding{UTF-8}
\title{
Functions to return fit statistics
}
\description{
\code{fitstats} produces sum of squared errors, Chi-squared statistic and Freeman-Tukey statistic used in parboot GOF tests throughout the book, starting with AHM1 - 7.5.4. \code{fitstats2} produces the above, plus corresponding statistics based on total detection frequency across replicate surveys; see AHM1 - 7.9.3.
}
\usage{
fitstats(fm)

fitstats2(fm)
}

\arguments{
  \item{fm}{
A fitted model object as returned by \pkg{unmarked}.
}
}

\value{
For \code{fitstats}, a named vector of length 3 with sum of squared errors (\code{SSE}), Chi-squared statistic (\code{Chisq}) and Freeman-Tukey (\code{freemanTukey}).

For \code{fitstats2}, a named vector of length 6 with the same information plus corresponding statistics based on numbers(\code{SSE.n}), (\code{Chisq.n}) and (\code{freemanTukey.n}).

}
\references{
Kéry & Royle (2016) \emph{Applied Hierarchical Modeling in Ecology} AHM1 - 7.5.4 and 7.9.3.
}
\author{
Marc Kéry & Andy Royle
}
