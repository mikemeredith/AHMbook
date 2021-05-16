
# AHM2 Section 7.2.2

# Named fp.modSel in earlier drafts.

# modSelFP <- function(mod.list){   # Thanks to Chris Sutherland!
  # (the same of the 'one-s**t hypothesis')
  # Model <- names(mod.list)
  # nPars <- sapply(mod.list, function(x)length(x@opt$par))
  # nll <- sapply(mod.list, function(x)x@opt$value)
  # AIC <- round(sapply(mod.list, function(x)x@AIC),2)
  # dAIC <- round(AIC - min(AIC), 2)
  # AICwt <- round(exp(-0.5 * dAIC)/sum(exp(-0.5 * dAIC)), 2)
  # modTab <- data.frame(nPars, AIC, dAIC, AICwt, row.names = Model)[order(AIC),]
  # modTab$cuWt <- cumsum(modTab$AICwt)  # Do this after sorting.
  # return(modTab)
# }

modSelFP <- function(mod.list){
  message("Please use the functions in 'unmarked' for this:\n  modSel(fitList(fits=mod.list))")
  unmarked::modSel(unmarked::fitList(fits=mod.list))
}

