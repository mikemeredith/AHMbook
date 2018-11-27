\name{graphSSM}
\alias{graphSSM}
\encoding{UTF-8}

% former name 'graph.ssm2'

\title{
Plot the results for multivariate time series of counts
}
\description{
Plot the observed time-series and the estimated population trajectories for state-space models (SSMs).
}
\usage{
graphSSM(ssm, C)
}
\arguments{
  \item{ssm}{
An object of class \code{jagsUI} with the results of fitting the SSM model.
}
  \item{C}{
The original count data: an individuals x year matrix of counts, no missing values.
}
}
\value{
No output, the function is called for plotting.
}
\references{
Kéry & Royle (2019) \emph{Applied Hierarachical Modeling in Ecology} Section 12.6

Kéry and Schaub (2012) \emph{Bayesian population analysis using WinBUGS - a hierarchical perspective}, Academic Press.
}
\author{
Adapted for multivariate SSMs from the SSM graphing function \code{graph.ssm} in chapter 5 in Kéry & Schaub (2012).
}

\examples{
# Would need the fitted model object to do this.
}
