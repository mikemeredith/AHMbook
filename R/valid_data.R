# AHMbook 18.5bis.2 p.52-55 in MS "Chapter_18_FINAL.docx".

# The function 'valid.data' appears in the supplementary materials (Appendix 2) of
# Chambert, T., Waddle, J.H., Miller, D.A.W., Walls, S.C., & Nichols, J.D. (2017) A new framework for analysing automated acoustic species-detection data: occupancy estimation and optimization of recordings post-processing. Methods in Ecology and Evolution, 9, 560-570.

# 'valid_data' performs the same function, but with new code by Mike Meredith, 2019-03-22.

# Both implementations select detections to validate in a series of rounds until the desired number, n.valid, has been reached. In the last round, not all candidate sites can be included without exceeding n.valid; valid.data takes the first sites in the list, valid_data takes a random sample.

### Arguments:
# N = vector of ALL detection counts (TP + FP)
# tp = vector of true positive (TP) detection counts
# n.valid = NUMBER of detections to be validated (if prop.valid=FALSE)
# n.valid = PROPORTION of detections to be validated (if prop.valid=TRUE)

# Returns a list with 2 components:
# n : the number of detections validated at each site
# k : the number of detections checked and found to be valid at each site

valid_data <- function(N, tp, n.valid, prop.valid=FALSE) {

  if(prop.valid){
    n.valid <- round(n.valid*sum(N))
  }
# -------------- Check and fix input -----------------------
  stopifnotEqualLength(N, tp)
  stopifnot(all(N >= tp))
  n.valid <- round(n.valid[1])
  stopifNegative(n.valid)
  # ------------------------------------------------------------
  if(n.valid > sum(N))
    warning("n.valid is greater than the total number of detections\n",
      "ALL will be validated", call.=FALSE)
  if(n.valid >= sum(N))
    return(list(n=N, k=tp))

  nsites <- length(N)

  # We will treat detections 1:tp[i] as true, the rest false.
  # We will randomise the selection of detections to validate.
  # (We do not need to randomise both.)

  # We will need to draw > 1 value per site, without drawing any twice,
  #   so we decide now on the (random) order of the draws for each site.
  # Each round of validation corresponds to one column of 'order'.
  order <- matrix(0, nrow=nsites, ncol=max(N))
  for(i in 1:nsites) {
    if(N[i] > 0)
      order[i, 1:N[i]] <- sample.int(N[i])
  }
  stopifnot(all(rowSums(order > 0) == N)) # check

  # Do successive rounds of validation until we have enough:
  n <- k <- numeric(nsites)  # vectors of zeros
  wanted <- n.valid - sum(n) # how many do we need to check?
  for(i in 1:max(N)) {
    todo <- which(order[, i] > 0)  # which sites have i (or more) detections
    if(length(todo) > wanted)      # if more than we want, subsample
      todo <- sample(todo, size=wanted)
    n[todo] <- n[todo] + 1
    k[todo] <- k[todo] + (order[todo,i] <= tp[todo])  # these are the good detections
    wanted <- n.valid - sum(n)   # how many do we still need to check?
    if(wanted == 0)
      break
  }

  return(list(n=n, k=k))
}
