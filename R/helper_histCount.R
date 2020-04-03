
# A helper function to plot a 'histogram' of counts (non-negative integers).

# Used in simNmix (2 places), data.fn, sim.fn, simComm.

# NOT EXPORTED

# Mike Meredith, 2016-12-25, updated 2020-04-03

histCount <- function(C1, C2 = NULL, x0 = FALSE, nbmax = 30, # max number of bins
    color = c('skyblue', adjustcolor('red', 0.5)), border='white', bty='n',
    ylab="Frequency", xlab="Count", main="", ...) {
# Plots a 'histogram' of 1 or 2 sets of counts (non-negative integers).
#
# C1 : a vector of non-negative integers; negative values silently ignored.
# C2 : an optional second vector of non-negative integers, or NULL.
# x0 : if TRUE, the x axis begins at zero, otherwise the lowest frequency.
# color : a length 2 vector of colours for the two histograms; the second colour
#   should be semi-transparent; can be scalar if C2 = NULL.
# border : colour of the border around bars, or NA for no border.
# other arguments as usual for plotting functions.

  rng <- range(C1, C2)
  if(x0)
    rng[1] <- 0
  range <- diff(rng) + 1
  if(range <= nbmax) {
    bwidth=1
  } else {
    bwidth <- range %/% nbmax + 1
  }
  nb <- ceiling(range/bwidth)
  br <- seq(min(C1,C2), max(C1,C2)+bwidth, by=bwidth) - 0.5

  H1 <- hist(C1, breaks = br, plot=FALSE)
  ymax <- max(H1$counts)
  if(!is.null(C2)) {
    H2 <- hist(C2, breaks = br, plot=FALSE)
    ymax <- max(ymax, H2$counts)
  }
      
  plot(H1, col=color[1], border=border,yaxs='i', ylim=c(0, ymax),
      xlim=rng, bty=bty, ylab=ylab, xlab=xlab, main=main)
  if(!is.null(C2)) {
    plot(H2, col=color[2], border=border,add=TRUE)
    plot(H1, breaks = br, col=NA, border=border,add=TRUE)  # replot borders of first histogram
  }
  axis(2)  # replot axis
  segments(rng[1], 0, rng[2], 0)  # do line at foot of bars
}
