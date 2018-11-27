
# Centre and scale a vector or array and return an object of the same class.
# For an array, the mean and SD of the whole array is used.

# Code modified from the function base::scale.

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
  if (is.numeric(center))
    attr(x, "scaled:center") <- center
  if (is.numeric(scale))
    attr(x, "scaled:scale") <- scale
  x
}

