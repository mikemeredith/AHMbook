\encoding{UTF-8}
\name{valid_data}
\alias{valid_data}
\title{
Validates simulated observational data against data-generating values.}
\description{
Where species identifications (from photos, recorded calls, etc) is doubtful, it is sometimes feasible for a subset of the data to be validated by experts. Chambert et al (2017) used simulations to investigate the effectiveness of such partial validation. Where a data set with false positives has been simulated and the true status is known, \code{valid_data} compares a subset with truth and returns the number of detections checked per site and the number checked and found to be valid. The function attempts to validate equal numbers of detections at each site (not equal proportions) subject to the actual number of detections and the total number to be validated; if the number of detections is small, all will be validated, otherwise a random sample is used; see Examples.

To recreate the data sets used in the book with R 3.6.0 or later, include \code{sample.kind="Rounding"} in the call to \code{set.seed}. This should only be used for reproduction of old results.
}
\usage{
valid_data(N, tp, n.valid, prop.valid=FALSE)

}
\arguments{
  \item{N}{
a vector with the number of detections for each site, including true positives and false positives.
}
  \item{tp}{
a vector with the number of true positive detections for each site.
}
  \item{n.valid}{
the number of detections to be validated (if \code{prop.valid=FALSE)} or the proportion of detections to be validated (if \code{prop.valid=TRUE)}.}
  \item{prop.valid}{
logical; determines how \code{n.valid} is interpreted.
}
}
\value{
A list with components:

\item{n }{a vector with the number of detections checked for each site.}

\item{k }{a vector with the number of detections checked and found to be valid (true positives) for each site.}
}
\references{
Chambert, T., Waddle, J.H., Miller, D.A.W., Walls, S.C., & Nichols, J.D. (2017) A new framework for analysing automated acoustic species-detection data: occupancy estimation and optimization of recordings post-processing. \emph{Methods in Ecology and Evolution}, 9, 560-570.
%% "analysing" is correct, this is a British journal.

Kéry, M. & Royle, J.A. (2021) \emph{Applied Hierarchical Modeling in Ecology} AHM2 - 7.6.2.
}

\author{
Mike Meredith.
}
\examples{
# Generate some data for 100 sites
set.seed(42, kind="Mersenne-Twister")
z <- rbinom(100, 1, 0.7)  # z[i] == 1 if species present at site 1, 0 otherwise
tp <- rpois(100, 3*z)     # Number of true detections, 0 if species not present
N <- tp + rpois(100, 0.5) # Add false positives

# Validate a subset of 150 detections
out <- valid_data(N, tp, 150)
head(tmp <- cbind(z=z, tp=tp, N=N, n=out$n, k=out$k), 10)
colSums(tmp)

# Plot the number validated vs all detections:
graphics::sunflowerplot(N, out$n, xlab="Number of detections at each site",
  ylab="Number validated")
# For sites with <= 2 detections, all are validated; otherwise at least 2 are validated,
#   with a 3rd draw from randomly selected sites to bring the total to 150.
}
