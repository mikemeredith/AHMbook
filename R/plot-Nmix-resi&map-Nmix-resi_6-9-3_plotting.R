# Functions for the book Applied Hierarchical Modeling in Ecology (AHM)
# Marc Kery & Andy Royle, Academic Press, 2016.

# plot.Nmix.resi and map.Nmix.resi - AHM1 section 6.9.3 p261


# Function to produce some residual plots in AHM1 Section 6.9.3
plot_Nmix_resi <- function(fmP, fmNB, fmZIP){
# Function does diagnostic plots for one Nmix model fitted with all three
#   mixture distributions currently availabe in unmarked:
#   Poisson, negative binomial and zero-inflated Poisson
# For each, fitted values vs. observed data and
#   residuals vs. fitted values are plotted.

# Plot fitted vs. observed data
op <- par(mfrow = c(2,3), mar = c(4,4,2,2), cex = 1.2) ; on.exit(par(op))
tmp1 <- range(c(fitted(fmP), fitted(fmNB), fitted(fmZIP)), na.rm = T)
limits1 = round(c(tmp1[1], tmp1[2]))
tmp2 <- range(c(residuals(fmP), residuals(fmNB), residuals(fmZIP)), na.rm = T)
limits2 = round(c(tmp2[1], tmp2[2]))

plot(fitted(fmP)~ fmP@data@y, xlab = "Observed data", ylab = "Fitted values (P)", frame = FALSE, ylim = limits1)
abline(0,1, lwd = 3 )
abline(lm(c(fitted(fmP))~ c(fmP@data@y)), col = "blue", lwd = 3)
plot(fitted(fmNB)~ fmP@data@y, xlab = "Observed data", ylab = "Fitted values (NB)", frame = FALSE, ylim = limits1)
abline(0,1, lwd = 3)
abline(lm(c(fitted(fmNB))~ c(fmP@data@y)), col = "blue", lwd = 3)
plot(fitted(fmZIP)~ fmP@data@y, xlab = "Observed data", ylab = "Fitted values (ZIP)", frame = FALSE, ylim = limits1)
abline(0,1, lwd = 3)
abline(lm(c(fitted(fmZIP)) ~ c(fmP@data@y)), col = "blue", lwd = 3)

# Plot residuals vs. fitted values
plot(residuals(fmP)~ fitted(fmP), xlab = "Fitted values (P)", ylab = "Residuals", frame = FALSE, xlim = limits1, ylim = limits2)
abline(h = 0, lwd = 2)
abline(lm(c(residuals(fmP)) ~ c(fitted(fmP))), col = "blue", lwd = 3)
plot(residuals(fmNB)~ fitted(fmNB), xlab = "Fitted values (NB)", ylab = "Residuals", frame = FALSE, xlim = limits1, ylim = limits2)
abline(h = 0, lwd = 2)
abline(lm(c(residuals(fmNB)) ~ c(fitted(fmNB))), col = "blue", lwd = 3)
plot(residuals(fmZIP)~ fitted(fmZIP), xlab = "Fitted values (ZIP)", ylab = "Residuals", frame = FALSE, xlim = limits1, ylim = limits2)
abline(h = 0, lwd = 2)
abline(lm(c(residuals(fmZIP)) ~ c(fitted(fmZIP))), col = "blue", lwd = 3)
}

# ..................................................................................................................



# Function to produce a map of the residuals in AHM1 Section 6.9.3
map.Nmix.resi <- function(fm, x, y){
# Function produces a map of the mean residuals from an N-mixture model
#    object named fm, which was fit by function pcount in unmarked
# Function arguments are the fitted model object and the x and y coordinates
#    of every site
mean.resi <- apply(residuals(fm), 1, mean, na.rm = TRUE)
mean.resi[mean.resi == "NaN"] <- mean(mean.resi, na.rm = TRUE)
spdata <- data.frame(residuals = mean.resi, x = x, y = y)
sp::coordinates(spdata) <- c("x", "y")
plot(sp::bubble(spdata, "residuals", col = c("blue", "red"), main = paste("Average residuals of fitted N-mixture model")))
}

