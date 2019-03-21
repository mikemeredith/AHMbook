\encoding{UTF-8}
\name{valid_data}
\alias{valid_data}
\title{
Validates simulated observational data against data-generating values.}
\description{
In cases where species identifications (from photos, recorded calls, etc) is doubtful, it is feasible for a subset of the data to be validated by experts. Chambert et al (2017) used simulations to investigate the effectiveness of such partial validation. Where a data set with false positives has been simulated and the true status is known, \code{valid_data} compares a randomly-selected subset with truth and returns the result.% Not sure what it returns right now!
}
\usage{
valid_data(N, tp, n.valid, prop.valid=FALSE)

}
\arguments{
  \item{N}{
data matrix of ALL detection counts (true positives and false positives).
}
  \item{tp}{
data matrix of true positive (TP) detection counts.
}
  \item{n.valid}{
the number of detections to be validated (if \code{prop.valid=FALSE)} or the proportion of detections to be validated (if \code{prop.valid=TRUE)}.}
  \item{prop.valid}{
logical; determines how \code{n.valid} is interpreted.
}
}
\value{
A list with components n and k.% More later!!
}
\references{
Chambert, T., Waddle, J.H., Miller, D.A.W., Walls, S.C., & Nichols, J.D. (2017) A new framework for analysing automated acoustic species-detection data: occupancy estimation and optimization of recordings post-processing. \emph{Methods in Ecology and Evolution}, 9, 560-570.

Kéry & Royle (2019) \emph{Applied Hierarachical Modeling in Ecology} Section 18...
}

\author{
Thierry Chambert and David A. W. Miller. % David L. Miller is at St Andrews and does distance sampling things.
% Their appendix 1 (BUGS code) starts with a long legalistic header from USGS, should we include that somewhere?

Function included in Appendix 2 to Chambert et al (2017).
}
\examples{

## Maybe not necessary, would need a simulated data set.

}
