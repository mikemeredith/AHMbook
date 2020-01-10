\name{plot_Nmix_resi}
\alias{plot_Nmix_resi}
\alias{plot.Nmix.resi}
\encoding{UTF-8}
\title{
Produce some residual plots
}
\description{
Function does diagnostic plots for one binomial N-mixture model fitted with all three mixture distributions currently available in package \pkg{unmarked}: Poisson, negative binomial and zero-inflated Poisson. For each, fitted values vs. observed data and residuals vs. fitted values are plotted.
}
\usage{
plot_Nmix_resi(fmP, fmNB, fmZIP)
}
\arguments{
  \item{fmP}{
Fitted model object for Poisson distribution.
}
  \item{fmNB}{
Fitted model object for negative binomial distribution.
}
  \item{fmZIP}{
Fitted model object for zero-inflated Poisson distribution.
}
}

\value{
None, used for its plotting effect.
}
\references{
Kéry, M. & Royle, J.A. (2016) \emph{Applied Hierarchical Modeling in Ecology} AHM1 - 6.9.3.

}
\author{
 Marc Kéry & Andy Royle
}
\note{
This function appears in the text as \code{plot.Nmix.resi}; renamed here to avoid confusion with generic \code{plot} functions.
}
