
# AHM2 section 3.4.1

# Define a function to create an m-array based on capture-histories (CH)

# Modified from Kery & Schaub (2012), a couple of loops replaced
#   with vector operations by Mike.

ch2marray <- function(CH){
  CH <- as.matrix(CH)  # might be a data frame
  nind <- nrow(CH)
  n.occasions <- ncol(CH)
  m.array <- matrix(data = 0, ncol = n.occasions+1, nrow = n.occasions)
    # First column and last row will be removed later
    # Last col is for number never-seen-again

  # Calculate the number of released individuals at each time period
  m.array[,1] <- colSums(CH)
  for (i in 1:nind){
    pos <- which(CH[i,]!=0) # When was animal caught?
    for (z in seq_along(pos[-1])) { # Does nothing if length(pos) == 1
      m.array[pos[z], pos[z+1]] <- m.array[pos[z], pos[z+1]] + 1
    } #z
  } #i
  # Calculate the number of individuals that is never recaptured
  m.array[, n.occasions+1] <- m.array[, 1] - rowSums(m.array[, -1])
  # Remove last row (releases on last occasion will never be recaptured)
  #   and 1st col (no REcaptures on 1st occasion).
  out <- m.array[-n.occasions, -1]
  return(out)
}
