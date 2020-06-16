
# Function to deal with errors arising in plotting stuff with sim* functions.
# Plotting errors should not cause the 'sim*()' function to halt but should produce
#   a sensible warning.

# not exported

# 'tryError' is the output from a 'try' call wrapped around plotting code

tryPlotError <- function(tryError) {
  msg <- "Plotting of output failed"
  msg2 <- attr(tryError, "condition")$message
  if(!is.null(msg2)) {
    if(msg2 == "figure margins too large")
      msg2 <- "the plotting window is too small."
    msg <- paste(msg, msg2, sep="\n   ")
  }
  if(Sys.getenv("RSTUDIO") == "1")
    msg <- paste(msg, "Try calling 'dev.new()' before the 'sim*' function.", sep="\n   ")
  warning(msg, call. = FALSE)
}
