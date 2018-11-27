\encoding{UTF-8}
\name{standardize}
\alias{standardize}
\title{
Scaling and centering of vectors, matrices and arrays
}
\description{
Maps a numeric variable to a new object with the same dimensions. Typically used to standardise a covariate to mean 0 and SD 1.
}
\usage{
standardize(x, center = TRUE, scale = TRUE)
}
\arguments{
  \item{x}{
a numeric vector, matrix or multidimensional array; \code{NA} and \code{NaN} are allowed; \code{Inf} and \code{-Inf} will produce all-\code{NaN} output if either \code{center} or \code{scale} is \code{TRUE}.
}
  \item{center}{
either a logical or a numeric value of length 1.
}
  \item{scale}{
either a logical or a numeric value of length 1.
}
}
\details{
\code{standardize} differs from \code{\link{scale}} by (1) accepting multidimensional arrays but not data frames; (2) not standardizing column-wise but using a single value to center or to scale; (3) if \code{x} is a vector, the output will be a vector (not a 1-column matrix).

Centering is performed before scaling.

If \code{center} is numeric, that value will be subtracted from the whole object. If logical and TRUE, the mean of the object (after removing NAs) will be subtracted.

If \code{scale} is numeric, the whole object will be divided by that value. If logical and TRUE, the standard deviation of the object (after removing NAs) will be used; this may not make sense if \code{center = FALSE}.
}
\value{
A numeric object of the same dimensions as \code{x} with the standardized values. NAs in the input will be preserved in the output.

For the default arguments, the object returned will have mean approximately zero and SD 1. (The mean is not exactly zero as scaling is performed after centering.)
}

\author{
Mike Meredith, after looking at the code of \code{base::scale}.
}
\examples{
# Generate some fake elevation data:
elev <- runif(100, min=100, max=500)
mean(elev)
sd(elev)
e <- standardize(elev)
mean(e)
sd(e)

# Standardize so that e=0 corresponds to exactly 300m and +/- 1 to
#   a change of 100m:
e <- standardize(elev, center=300, scale=100)
mean(e)
mean(elev) - 300
range(e)
range(elev) - 300
}
