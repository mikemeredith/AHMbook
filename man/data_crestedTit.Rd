\name{crestedTit}
\alias{crestedTit}
\encoding{UTF-8}
\docType{data}
\title{
Data from the Swiss Breeding Bird Survey for crested tits from 1999 to 2015
}
\description{
The Swiss breeding bird survey ("Monitoring Häufige Brutvögel" MHB) has monitored the populations of 150 common species since 1999. The MHB sample consists of 267 1-km squares that are laid out as a grid across Switzerland. Fieldwork is conducted by about 200 skilled birdwatchers, most of them volunteers. Avian populations are monitored using a simplified territory mapping protocol, where each square is surveyed up to three times during the breeding season (only twice above the tree line). Surveys are conducted along a transect that does not change over the years.

The data frame \code{crestedTit} has the data for crested tits \emph{Parus cristatus} from 1999 to 2015. There are some missing values: see Details.
}

\usage{data("crestedTit")}

\format{
\code{crestedTit} is a data frame with 267 rows and 133 columns:

\describe{
  \item{specid }{a numeric species ID based on phylogeny.}
  \item{latid }{a 6-letter abbreviation of the Latin name.}
  \item{coordx, coordy }{the x and y coordinates of the centre of the quadrat; the coordinate reference system intentionally not specified.}
  \item{elev }{the mean elevation of the quadrat, m.}
  \item{forest }{percentage forest cover}
  \item{nsurveys }{the number of replicate surveys \emph{planned} in the quadrat; above the tree-line 2, otherwise 3.}
   \item{X1999 to X2016 }{total number of territories observed for the year.}
   \item{date991 to date163 }{Julian date of the survey, 1 = 1st January; the first 2 digits indicate the year and the 3rd digit the survey.}
   \item{dur991 to dur163 }{duration of the survey, minutes; the first 2 digits indicate the year and the 3rd digit the survey.}
  }
}
\details{

Missing values in the \code{date} array indicate that the corresponding survey was not carried out.

On 26 occasions when surveys were carried out, the duration was not recorded, resulting in additional NAs in the \code{dur} array.

A new method for recording breeding territories was introduced in 2004, but the old protocol was in use at some sites until 2013. %Surveys with the old protocol have the counts shown as NA in the \code{count} array.
}

\source{
Swiss Ornithological Institute
}
\references{
Kéry & Royle (2019) \emph{Applied Hierarachical Modeling in Ecology} chapter 12.
}
\examples{
data(crestedTit)
str(crestedTit)


}
\keyword{datasets}
