# Functions for the book Applied Hierarchical Modeling in Ecology (AHM)
# Marc Kéry & Andy Royle, Academic Press, 2016.

# ppc.plot - section 68 p253


# Function to plot results from posterior predictive check in AHM section 6-8,
# for a fitted model object with JAGS, as in that section

ppc.plot <- function(fm){
par(mfrow = c(2,2), mar = c(5,5,3,2), cex.lab = 1.3, cex.axis = 1.3)
# Function plots results from posterior predictive check
#   in AHM section 6-8 for a fitted model object with JAGS
fit.a <- fm$sims.list$fit.actual    # Extract posterior samples
fit.s <- fm$sims.list$fit.sim
ch <- fm$sims.list$c.hat
lims <- c(min(c(fit.a, fit.s)), max(c(fit.a, fit.s)))
hist(fit.a, breaks = 100, col = "grey", main = "", xlab = "Fit statistic actual data")
hist(fit.s, breaks = 100, col = "grey", main = "", xlab = "Fit statistic simulated data")
hist(ch, breaks = 100, col = "grey", main = "", xlab = "Lack of fit ratio (c-hat)")
title(paste("c-hat =", round(mean(ch),2)))
plot(fit.a[fit.a >= fit.s], fit.s[fit.a >= fit.s], xlab = "Fit statistic actual data",
ylab = "Fit statistic simulated data", col = "blue", xlim = lims, ylim = lims, frame = F)
title(paste("bpv (proportion red) =", round(mean(fit.s>fit.a),2)))
points(fit.a[fit.a < fit.s], fit.s[fit.a < fit.s], col = "red")
abline(0,1)
}

