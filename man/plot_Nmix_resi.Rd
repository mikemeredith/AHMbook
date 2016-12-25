\name{plot_Nmix_resi}
\alias{plot_Nmix_resi}
\encoding{UTF-8}
\title{
Produce some residual plots
}
\description{
Function does diagnostic plots for one Nmix model fitted with all three mixture distributions currently availabe in package \pkg{unmarked}: Poisson, negative binomial and zero-inflated Poisson. For each, fitted values vs. observed data and residuals vs. fitted values are plotted.
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
Kéry & Royle (2016) \emph{Applied Hierarachical Modeling in Ecology} Section 6.9.3

}
\author{
 Marc Kéry
}
\note{
This function appears in the text as 'plot.Nmix.resi'; renamed to avoid confusion with generic 'plot' functions.
}

\examples{
##
}
