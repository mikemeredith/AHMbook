
# Functions 'standardize' and 'standardize2match'

# Centre and scale a vector or array and return an object of the same class.
# For an array, the mean and SD of the whole array is used.

standardize <- function (x, center = TRUE, scale = TRUE) {
  if (!is.numeric(x))
    stop("'x' must be a numeric vector or array.", call. = FALSE)
  if (length(center) != 1)
      stop("'center' must be logical or numeric of length 1.", call. = FALSE)
  if (length(scale) != 1)
      stop("'scale' must be logical or numeric of length 1.", call. = FALSE)

  if (is.logical(center)) {
    if (center) {
      center <- mean(x, na.rm = TRUE)
      x <- x - center
    }
  } else {
    if (!is.numeric(center))
       stop("'centre' must be numeric or logical.", call. = FALSE)
    x <- x - center
  }
  if (is.logical(scale)) {
    if (scale) {
      scale <- sd(x, na.rm=TRUE)
      x <- x / scale
    }
  } else {
    if (!is.numeric(scale))
       stop("'scale' must be numeric or logical.", call. = FALSE)
    x <- x / scale
  }
  return(x)
}
#........................................................................

# Standardize a new numeric object to the same mean and sd as 
#   existing output from 'standardize'
standardize2match <- function (x, y) {
  if (!is.numeric(x) || !is.numeric(x))
    stop("'x' and 'y' must be a numeric vectors or arrays.", call. = FALSE)
  return((x - mean(y, na.rm=TRUE)) / sd(y, na.rm=TRUE))
}

