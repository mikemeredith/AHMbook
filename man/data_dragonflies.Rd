\name{dragonflies}
\alias{dragonflies}
\alias{pop}
\alias{sex}
\alias{wing}
\alias{body}
\alias{mites}
\alias{color}
\alias{damage}
\docType{data}
\encoding{UTF-8}
\title{
Fictional data for dragonflies
}
\description{
A toy data set with fictional data for 9 individuals of the dragonfly \emph{Onychogomphus uncatus} from 3 populations in the Spanish Pyrenees.
}
\usage{
data(dragonflies)}

\format{
  The format is seven vectors of length 9:
\describe{
  \item{pop}{a factor indicating which population the individual was drawn from.}
  \item{sex}{a factor indicating the sex of each individual.}
  \item{wing}{wingspan.}
  \item{body}{body length.}
  \item{mites}{number of mites (ectoparasites) counted on each individual.}
  \item{color}{proportion of body yellow as opposed to black.}
  \item{damage}{the number of wings (out of 4) damaged.}
}
}
\source{
Fictitious data.
}
\references{
Kéry, M. & Royle, J.A. (2016) \emph{Applied Hierarchical Modeling in Ecology} AHM1 - 3.1
}
\examples{
data(dragonflies)
lm(wing ~ pop + body)

}
\keyword{datasets}
