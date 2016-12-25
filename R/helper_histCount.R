
# A helper function to plot a 'histogram' of counts (non-negative integers).

# Used in simNmix (2 places), data.fn, sim.fn, simComm.

# NOT EXPORTED

# Mike Meredith, 2016-12-25

histCount <- function(C1, C2 = NULL, x0 = FALSE,
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

  CC <- list(C1, C2)
  freq <- vector('list', 2)
  xx <- vector('list', 2)
  for(i in 1:2) {
    if(is.null(CC[[i]]))
      next
    tab <- tabulate(round(CC[[i]]) + 1) # +1 because tabulate ignores zero
    xx[[i]] <- seq(min(which(tab > 0)), length(tab), 1) - 1
    freq[[i]] <- tab[xx[[i]] + 1]
  }
  xlim <- range(xx) + c(-0.5, +0.5)
  if(x0)
    xlim[1] <- -0.5
  ylim <- range(0, freq)*1.04 # add 'head room' to ylim, as yaxs='i'.

  plot(xlim, ylim, type='n', yaxs='i',
    bty=bty, ylab=ylab, xlab=xlab, main=main, ...)
  for(i in 1:2) {
    if(is.null(xx[[i]]))
      next
    rect(xx[[i]]-0.5, 0, xx[[i]]+0.5, freq[[i]], col=color[i],
      border=border, ...)
  }
  if(!is.na(border) && !is.null(xx[[2]])) # replot borders of 1st histogram
    rect(xx[[1]]-0.5, 0, xx[[1]]+0.5, freq[[1]], col=NA, border=border, ...)
  segments(xlim[1], 0, xlim[2], 0)  # do line at foot of bars
}

# Some examples
if(FALSE) {
  C1 <- rpois(100, 8)
  C2 <- rpois(100, 0.8)
  histCount(C1, C2)
  histCount(C1, C2, bty="l", las=1, lwd=5, main="Cool histogram")
  histCount(C1,  , las=1, main="Cool histogram")
  histCount(C1 + 10,  , las=1, main="Cool histogram", x0=TRUE)
}

