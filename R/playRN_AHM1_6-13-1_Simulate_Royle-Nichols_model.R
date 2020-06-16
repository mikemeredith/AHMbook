# Functions for the book Applied Hierarchical Modeling in Ecology (AHM)
# Marc Kery & Andy Royle, Academic Press, 2016.

# playRN - AHM1 section 6.13.1 p301

# Function to play Royle-Nichols model
#   (introduced in AHM1 Section 6.13.1)

playRN <- function(M = 267, J = 3, mean.abundance = 1, mean.detection = 0.3,
    show.plots = TRUE, verbose = TRUE){
# Function generates replicated count data under the Nmix model of Royle (2004),
#   then 'degrades' the data to detection/nondetection and fits the RN model
#   (Royle & Nichols 2003) using unmarked and estimates site-specific abundance.
#   Requires function simNmix and package unmarked.
#
# devAskNewPage(ask = FALSE) ## leave it as it is!
#
# Simulate Nmix data under a range of abundance levels
data <- simNmix(nsites = M, nvisits = J, mean.lam = mean.abundance, mean.p = mean.detection,
  beta2.lam = 1, beta3.p = -1, beta.p.survey = -1, show.plots = FALSE, verbose = verbose)
# Turn counts into detection/nondetection data
y <- data$C          # Copy counts C into y
y[y>0] <- 1          # Turn counts >0 into 1
# Load unmarked, format data and summarize
umf <- unmarkedFrameOccu(y=y, siteCovs= data.frame(cov2 = data$site.cov[,2],
  cov3 = data$site.cov[,3]), obsCovs = list(obscov = data$survey.cov))
# Fit data-generating model
fm <- occuRN(~cov3+obscov ~cov2, data=umf)
# Estimate local abundance N and plot against true N (known in simulation)
Nest <- bup(ranef(fm, K = ), "mean")
if(show.plots) {
  # par(mfrow = c(1,1)) ## leave it alone
  tryPlot <- try( {
    plot(data$N, Nest, xlab = "True local abundance",
      ylab = "Estimated local abundance", frame = FALSE)
    abline(0,1, lwd = 3)                              # 1:1 line
    abline(lm(Nest ~ data$N), col = "blue", lwd = 3)  # Regression
  }, silent = TRUE)
  if(inherits(tryPlot, "try-error"))
    tryPlotError(tryPlot)
}
slope <- coef(lm(Nest ~ data$N))[2]               # Is 1 if model perfect
return(list(nsites = M, nvisits = J, coef = coef(fm), slope = slope))
}

