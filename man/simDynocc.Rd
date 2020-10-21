\name{simDynocc}
\alias{simDynocc}
\encoding{UTF-8}
\title{
Simulate detection/nondetection data under a wide variety of non-spatial dynamic occupancy models
}
\description{
Function to simulate detection/nondetection data under a general dynamic site-occupancy model, including:

* Annual variation in the probabilities of patch persistence, colonization and detection is specified by the bounds of a uniform distribution.

* One covariate is allowed to affect each parameter: a site covariate for psi1, site-by-year covariates for phi and gamma, and an observational covariate for p. Covariates are generated internally from uniform(-2, 2) distributions.

* Additional heterogeneity among sites in persistence and colonization or both.

* Additional detection heterogeneity at the site-, the survey, or the site-by-survey level, with the possibility of a temporal trend in this heterogeneity over the years. E.g., an annual trend in detection heterogeneity at the site or the survey level is specified by the first and second value, which correspond to the heterogeneity in the first and the last year. Hence, trend.sd.site = c(0, 1) will result in a linear trend in the magnitude of site heterogeneity in detection from 0 in the first year to 1 in the last year.

* Additional detection heterogeneity that varies over the survey (= occasion) according to a quadratic effect of occasion number (to model the typical phenology of an insect species for instance).

* Simulation of data under a BACI (before-after-control-impact) design, where some event happens in a given year and \emph{reduces} phi or gamma by a stated percentage (only reductions, no increases allowed!)

}
\usage{
simDynocc(nsites = 250, nyears = 10, nsurveys = 3, year.of.impact = NA,
  mean.psi1 = 0.4, beta.Xpsi1 = 0,
  range.phi = c(0.5, 1), sd.lphi.site = 0, impact.phi = 0, beta.Xphi = 0,
  range.gamma = c(0, 0.5), sd.lgamma.site = 0, impact.gamma = 0, beta.Xgamma = 0,
  sd.lphi.lgamma.site = 0,
  range.p = c(0.1, 0.9), beta.Xp = 0,
  range.beta1.survey = c(0, 0), range.beta2.survey = c(0, 0),
  trend.sd.site = c(0, 0), trend.sd.survey = c(0, 0),
  trend.sd.site.survey = c(0, 0), show.plots = TRUE)
}
\arguments{
  \item{nsites}{
Number of sites.
}
  \item{nyears}{
Number of years (or 'seasons').
}
  \item{nsurveys}{
Number of replicate surveys (= occasions) within a year.
}
  \item{year.of.impact}{
Year in which an impact happens that can affect phi and gamma (for BACI design), NA if no impact occurs; for the BACI design, \code{nyears} must be greater than 2, and year of impact cannot be the first or last year.
}
  \item{mean.psi1}{
average occupancy probability in first year.
}
  \item{beta.Xpsi1}{
coefficient of environmental covariate in probability of initial occupancy.
}
  \item{range.phi}{
bounds of uniform distribution from which annual probability of persistence is randomly drawn.
}
  \item{sd.lphi.site}{  % new 2019-03-01
SD of random site effect on persistence on the logit scale drawn from a normal distribution with mean zero.
}
  \item{impact.phi}{
negative effect in percent on annual phi (e.g., impact.phi = 20 means a 20\% reduction in phi); ignored if \code{year.of.impact = NA}.
}
  \item{beta.Xphi}{
coefficients of environmental covariate in probability of persistence.
}
  \item{range.gamma}{
bounds of uniform distribution from which annual probability of colonization is randomly drawn.
}
  \item{sd.lgamma.site}{  % new 2019-03-01
SD of random site effect on colonization on the logit scale drawn from a normal distribution with mean zero.
}
  \item{impact.gamma}{
negative effect in percent on annual gamma (e.g., impact.gamma = 20 means a 20\% reduction in gamma); ignored if \code{year.of.impact = NA}.
}
  \item{beta.Xgamma}{
coefficient of environmental covariate in probability of colonization.
}
  \item{sd.lphi.lgamma.site}{  % new 2019-03-01
SD of random site effect on persistence AND colonization on the logit scale drawn from a normal distribution with mean zero.
}
  \item{range.p}{
bounds of uniform distribution from which annual probability of detection is randomly drawn.
}
  \item{beta.Xp}{
coefficients of environmental covariate in probability of detection.
}
 \item{range.beta1.survey}{ %%%%
bounds of the uniform distribution from which the annual variation in the linear effect of the survey occasion (i.e., of survey 1-12 with \code{nsurveys = 12}) on detection (= product of availability and perceptibility) is randomly drawn.
}
  \item{range.beta2.survey}{ %%%%
the same for the quadratic effect of survey occasion.
}
  \item{trend.sd.site}{
initial and final values of sd of normal distribution to model logit-normal noise in p at the site level; a linear trend is assumed over time; if the two values are the same, a constant value is assumed.
}
  \item{trend.sd.survey}{
initial and final values of sd of normal distribution to model logit-normal noise in p \emph{only} at the 'survey' level; if they are different, a linear trend is assumed over time.
}
  \item{trend.sd.site.survey}{
initial and final values of sd of normal distribution to model logit-normal noise in p at the site/year/survey = ‘survey’ level; if they are different, a linear trend is assumed over time.
}
  \item{show.plots}{
If TRUE, plots of results are displayed; set to FALSE if running simulations.
}
}

\value{
A list with the values of the arguments input and the following additional elements:
  \item{impact }{a 0/1 vector of length (nyears-1) indicating if an impact applies to the interval}
  \item{BACI.effect.phi }{reduction in persistence due to impact, a vector of length nyears - 1}
  \item{BACI.effect.gamma }{reduction in colonization due to impact, a vector of length nyears - 1}
  \item{beta1 }{linear effect of occasion on the product of availability and detection, a vector of length nyears}
  \item{beta2 }{quadratic effect of occasion on the product of availability and detection, a vector of length nyears}
  \item{mean.phi }{mean persistence for each interval before application of any BACI effect, a vector of length nyears - 1}
  \item{mean.gamma }{mean colonization for each interval before application of any BACI effect, a vector of length nyears - 1}
  \item{mean.p }{mean detection probability for each year, a vector of length nyears}
  \item{psi }{annual occupancy for each site, a nsites x nyears matrix}
  \item{mean.psi }{average occupancy over sites, a vector of length nyears}
  \item{n.occ }{number of occupied sites, a vector of length nyears}
  \item{psi.fs }{finite-sample occupancy proportion, a vector of length nyears}
  \item{psi.app }{apparent occupancy over sites, a vector of length nyears}
  \item{z }{true occurrence state, a nsites x nyears matrix of 0/1}
  \item{phi }{persistence, a nsites x nyears-1 matrix}
  \item{gamma }{colonization, a nsites x nyears-1 matrix}
  \item{p }{detection probability, a nsites x nsurveys x nyears array}
  \item{y }{the observed detection history, a nsites x nsurveys x nyears array}
  \item{Xpsi1 }{covariate affecting initial occupancy, a vector of length nsites}
  \item{Xphi }{covariate affecting persistence, a nsites x nyears matrix}
  \item{Xgamma }{covariate affecting colonization, a nsites x nyears matrix}
  \item{Xp }{covariate affecting probability of detection, a nsites x nsurveys x nyears array}
  \item{eps.lphi.site }{site random effects on persistence, a vector of length nsites}
  \item{eps.lgamma.site }{site random effects on colonization, a vector of length nsites}
  \item{eps.lphi.lgamma.site }{site random effects on persistence and colonization, a vector of length nsites}
  \item{eps1 }{site random effects on detection, a vector of length nsites}
  \item{eps2 }{the survey random effects on detection, a vector of length nsurveys}
  \item{eps3 }{the site/survey/year random effects on detection, a nsites x nsurveys x nyears array}
}
\references{
Kéry, M. & Royle, J.A. (2021) \emph{Applied Hierarchical Modeling in Ecology} AHM2 - 4.

}
\author{
Marc Kéry & Andy Royle
}

\examples{
# Generate data with the default arguments and look at the structure
tmp <- simDynocc()
str(tmp)

# no annual variation in the parameters
str(data <- simDynocc(nsites = 250, nsurveys = 3, nyears = 10, mean.psi1 = 0.6,
    range.phi = c(0.7, 0.7), range.gamma = c(0.3, 0.3), range.p = c(0.5, 0.5)))
# a fully time-dependent model (with p constant within each primary period)
str(data <- simDynocc(mean.psi1 = 0.6, range.phi = c(0.5, 0.8),
   range.gamma = c(0.1, 0.5), range.p = c(0.1, 0.9)) )
# a time-constant model with four different covariates affecting the four parameters
str(data <- simDynocc(mean.psi1 = 0.6, beta.Xpsi1 = 1,
  range.phi = c(0.6, 0.6), beta.Xphi = 2, range.gamma = c(0.3, 0.3),
  beta.Xgamma = 2, range.p = c(0.2, 0.2), beta.Xp = -2) )
# seasonal variation in detection probability
str(data <- simDynocc(nsurveys = 12, mean.psi1 = 0.6,
  range.phi = c(0.6, 0.6), range.gamma = c(0.3, 0.3),
  range.p = c(0.5, 0.5), range.beta1.survey = c(-0.3, 0.4),
  range.beta2.survey = c(0, -0.7)) )
# now both yearly variation and effects of all covariates (including survey)
str( data <- simDynocc(mean.psi1 = 0.6, beta.Xpsi1 = 1,
   range.phi = c(0.6, 1), beta.Xphi = 2, range.gamma = c(0, 0.2),
   beta.Xgamma = 2, range.p = c(0.1, 0.9), beta.Xp = -2,
   range.beta1.survey = c(-0.4, 0.5), range.beta2.survey = c(0, -0.8)) )

# To add detection heterogeneity at the site level, you can do this:
str(data <- simDynocc(trend.sd.site = c(3, 3)) ) # No time trend
str(data <- simDynocc(trend.sd.site = c(1, 3)) ) # With time trend

# To add detection heterogeneity at the level of the survey, you can do this:
str(data <- simDynocc(trend.sd.survey = c(3, 3)) ) # No time trend
str(data <- simDynocc(trend.sd.survey = c(1, 3)) ) # With time trend

# To add detection heterogeneity at the level of the individual visit, you can do this:
str(data <- simDynocc(trend.sd.site.survey = c(3, 3)) ) # No trend
str(data <- simDynocc(trend.sd.site.survey = c(1, 3)) ) # With trend

# To simulate data under a BACI design, where an impact happens in year 10
str(data <- simDynocc(nsites = 250, nsurveys = 3, nyears = 20, year.of.impact = 10,
   impact.phi = 80, impact.gamma = 50) )

# And data where there is no detection error (i.e., with p = 1):
str( data <- simDynocc(range.p = c(1, 1)) )

}
