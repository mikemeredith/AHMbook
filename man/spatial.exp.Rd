\name{spatial.exp}
\alias{spatial.exp}
\encoding{UTF-8}

\title{
Generate random field with negative exponential correlation
}
\description{
Function creates Gaussian random field with negative exponential correlation and visualizes correlation and random field
}
\usage{
spatial.exp(variance = 1, theta = 1, size = 50, show.plot = TRUE)
}
\arguments{
  \item{variance}{
variance of field
}
  \item{theta}{
parameter governing spatial correlation (=1/phi) ("large theta means high correlation") Note that RMexp is specified in terms of phi = 1/theta
}
  \item{size}{
Number of pixels in either direction
}
  \item{show.plot}{
if TRUE, plots of the data will be displayed; set to FALSE if you are running many simulations or use inside other functions.
}
}
\value{
A  list with the following components:
 \item{variance }{the value of the variance argument}
 \item{theta }{the value of the theta argument}
 \item{field }{the random field variable, a size x size matrix}
 \item{grid }{the grid corresponding to field, equivalent to \code{expand.grid(1:size, 1:size)}}
}
\references{
Kéry & Royle (2016) \emph{Applied Hierarachical Modeling in Ecology} Chapter 20

}
\author{
Kéry & Royle
}

\examples{
# Generate data with the default arguments and look at the structure:
tmp <- spatial.exp()
str(tmp)
}
