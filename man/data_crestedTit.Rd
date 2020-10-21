\name{crestedTit}
\alias{crestedTit}
\encoding{UTF-8}
\docType{data}
\title{
Count data from the Swiss Breeding Bird Survey MHB for crested tits from 1999 to 2015
}
\description{
The Swiss breeding bird survey ("Monitoring Häufige Brutvögel" MHB) has monitored the populations of approximately 150 common species since 1999. The MHB sample consists of 267 1-km squares that are laid out as a grid across Switzerland. Fieldwork is conducted by about 200 skilled birdwatchers annually, most of them volunteers. Avian populations are monitored using a simplified territory mapping protocol (also called spot-mapping West of the Atlantic), where each square is surveyed up to three times during the breeding season (only twice above the tree line around 2000 m elevation). Surveys are conducted along a transect that does not change over the years and is 4-6 kms in length on average.

The data frame \code{crestedTit} has the data for Crested Tits \emph{Parus cristatus} from 1999 to 2015. There are some missing values: see Details.
}

\usage{data("crestedTit")}

\format{
\code{crestedTit} is a data frame with 267 rows and 131 columns:

\describe{
  \item{coordx, coordy }{the x and y coordinates of the center of the quadrat; the coordinate reference system intentionally not specified.}
  \item{elev }{the mean elevation of the quadrat, in meters.}
  \item{forest }{percentage forest cover}
  \item{nsurveys }{the number of replicate surveys \emph{planned} in the quadrat; above the tree-line 2, otherwise 3.}
   \item{1999 to 2016 }{total number of territories observed for the year.}
   \item{date991 to date163 }{Julian date of the survey, 1 = 1st January; the first 2 digits in the column name indicate the year and the 3rd digit the survey.}
   \item{dur991 to dur163 }{duration of the survey, in minutes; the first 2 digits in the column name indicate the year and the 3rd digit the survey.}
  }
}
\details{

Missing values in the \code{date} array indicate that the corresponding survey was not carried out.

On 26 occasions when surveys were carried out, the duration was not recorded, resulting in additional NAs in the \code{dur} array.

A new method for recording breeding territories was introduced in 2004, but the old protocol was in use at some sites until 2013. Surveys with the old protocol have the counts shown as NA in the \code{count} array.

See also Chapter 6 in Kéry & Royle (2016) for further description of the survey and the data it produces.
}

\source{
Swiss Ornithological Institute
}
\references{
Kéry & Royle (2016) \emph{Applied Hierarchical Modeling in Ecology} AHM1 - 6.

Kéry & Royle (2021) \emph{Applied Hierarchical Modeling in Ecology} AHM2 - 1.
}
\examples{
data(crestedTit)
str(crestedTit)

# Prepare objects needed for section AHM2 - 1.3 and following
C <- as.matrix(crestedTit[, 6:23]) # matrix of counts
year <- 1999:2016
datetmp <- as.matrix(crestedTit[, 24:77])  # matrix of Julian dates
datefull <- array(datetmp, c(nrow(datetmp), 3, ncol(datetmp)/3))
    # site x rep x year array
durtmp <- as.matrix(crestedTit[, 78:131])  # matrix of survey durations
durfull <- array(durtmp, c(nrow(durtmp), 3, ncol(durtmp)/3))
    # site x rep x year array

# Get mean date and duration of survey for each site and year
mdate <- apply(datefull, c(1,3), mean, na.rm=TRUE)
mdate[is.nan(mdate)] <- NA  # Replace NaN with NA
mdur <- apply(durfull, c(1,3), mean, na.rm=TRUE)
mdur[is.nan(mdur)] <- NA

# For Sec 1.5, we need standardised covariates with missing values imputed
elev.sc <- as.numeric(scale(crestedTit$elev))     # site elevation
forest.sc <- as.numeric(scale(crestedTit$forest)) # site forest cover

mean.date <- mean(mdate, na.rm=TRUE)         # mean survey date per site
sd.date <- sd(mdate, na.rm=TRUE)
mdate.sc <- (mdate - mean.date) / sd.date
mdate.sc[is.na(mdate.sc)] <- 0               # mean imputation

mean.dur <- mean(mdur, na.rm=TRUE)           # mean survey duration per site
sd.dur <- sd(mdur, na.rm=TRUE)
mdur.sc <- (mdur - mean.dur) / sd.dur
mdur.sc[is.na(mdur.sc)] <- 0                 # mean imputation

# For sections 1.6 and 1.7, we remove sites with no crested tits recorded,
#   or recorded in only one year:
nzero <- apply(C == 0, 1, sum, na.rm=TRUE) # number of zero years per site
sel <- nzero <= 1  # TRUE if site has 2 or more years of non-zero data

}
\keyword{datasets}
