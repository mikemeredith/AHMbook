# Functions for the book Applied Hierarchical Modeling in Ecology (AHM)
# Marc Kery & Andy Royle, Academic Press, 2016.

# HELPER FUNCTIONS

# e2dist

# Function to compute Euclidean distances
# (from the R package that goes along with the
#    Spatial Capture-Recapture book by Royle et al. (2014))

 e2dist <- function(x, y=NULL){
    if(ncol(x) != 2)
      stop("Argument 'x' must be a 2-column matrix or data frame.", call.=FALSE)
    if(is.null(y)) {
      y <- x
    } else {
      if(ncol(y) != 2)
        stop("Argument 'y' must be a 2-column matrix or data frame.", call.=FALSE)
    }
    i <- sort(rep(1:nrow(y), nrow(x)))
    dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
    matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = FALSE)
}
# ....................................................................

