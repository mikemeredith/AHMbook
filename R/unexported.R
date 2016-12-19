
# Unexported helper functions/utilities
# ==========---------------------------

# Plot 'histograms' of species richness.
# Used in simComm
# S1, S2 : vectors of numbers of species at each site; typically S1 = true, S2 = observed.
# colors : vector of colours for S1 and S2 respectively.
# offset : value added to S1 and subtracted for S2 so that the bars are distinguishable

plotSpeciesTables <- function(S1, S2, colors=c('red', 'blue'), lwd=2, xlim=range(S1, S2), xlab="", main="", ylab="Frequency", offset=0.2, ...) {
  S1t <- table(S1)
  S2t <- table(S2)
  plot(x=as.numeric(names(S1t))+offset, y=S1t, type='h', xlim=xlim, ylim=c(0, max(S1t,S2t)), lwd=lwd, col=colors[1],
    main=main, xlab=xlab, ylab=ylab, ...)
  axis(2, ...)
  points(as.numeric(names(S2t))-offset, S2t, type='h',lwd=lwd, col=colors[2], ...)
}
