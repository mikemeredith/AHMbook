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
Imaginary data for dragonflies
}
\description{
A toy data set with fictional data for 9 dragonflies \emph{Onychogomphus uncatus} from 3 populations in the Spanish Pyrenees.
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
  \item{mites}{number of mites (ectoparasites) counted.}
  \item{color}{proportion of body yellow as opposed to black.}
  \item{damage}{the number of wings (out of 4) damaged.}
}
}
\source{
Fictitious data.
}
\references{
Kéry & Royle (2016) \emph{Applied Hierarchical Modeling in Ecology} Section 3.1
}
\examples{
data(dragonflies)
lm(wing ~ pop + body)

}
\keyword{datasets}
