
# Sanity checks

# Not exported.

# These functions are very basic but simplify input checks while giving
#   informative error messages.
# Without checks, users may get mysterious error messages, eg,
#   "is.na() applied to non-(list or vector) of type 'closure'".
# base::stopifnot() is a start but error messages are still abstruse, eg,
#  "p >= 0 & p <= 1 is not TRUE".

stopifnotNumeric <- function(arg, allowNA=FALSE) {
  name <- deparse(substitute(arg))
  if(allowNA && all(is.na(arg))) {
    # do nothing
  } else {
    if(!allowNA && any(is.na(arg)))
      stop("Argument '", name, "' must not be NA or NaN.", call.=FALSE)
    if(!is.numeric(arg))
      stop("Argument '", name, "' must be numeric.", call.=FALSE)
  }
}

stopifnotInteger <- function(arg, allowNA=FALSE) {
  name <- deparse(substitute(arg))
  if(allowNA && all(is.na(arg))) {
    # do nothing
  } else {
    if(!allowNA && any(is.na(arg)))
      stop("Argument '", name, "' must not be NA or NaN.", call.=FALSE)
    if(!is.numeric(arg))
      stop("Argument '", name, "' must be numeric.", call.=FALSE)
    if(!all(arg%%1 == 0))
      stop("Argument '", name, "' must be integer (whole number).", call.=FALSE)
  }
}

stopifnotScalar <- function(arg, allowNA=FALSE) {
  name <- deparse(substitute(arg))
  if(length(arg) > 1)
    stop("Argument '", name, "' must be a single value.", call.=FALSE)
  if(allowNA && is.na(arg)) {
    # do nothing
  } else {
    if(!allowNA && is.na(arg))
      stop("Argument '", name, "' must not be NA or NaN.", call.=FALSE)
    if(!is.numeric(arg))
      stop("Argument '", name, "' must be numeric.", call.=FALSE)
  }
}

stopifnotLength <- function(arg, length) {
  name <- deparse(substitute(arg))
  if(length(arg) != length)
    stop("Argument '", name, "' must have length ", length, ".", call.=FALSE)
}

stopifnotProbability <- function(arg, allowNA=FALSE) {
  name <- deparse(substitute(arg))
  if(allowNA && all(is.na(arg))) {  # An all-NA vector is logical, but ok.
    # do nothing
  } else {
    if(!allowNA && any(is.na(arg)))
      stop("Argument '", name, "' must not be NA or NaN.", call.=FALSE)
    if(!is.numeric(arg))
      stop("Argument '", name, "' must be numeric.", call.=FALSE)
    if(any(arg < 0 | arg > 1, na.rm=TRUE))
      stop("Argument '", name, "' must be between 0 and 1.", call.=FALSE)
  }
}

stopifnotBetween <- function(arg, min, max, allowNA=FALSE) {
  name <- deparse(substitute(arg))
  if(allowNA && all(is.na(arg))) {  # An all-NA vector is logical, but ok.
    # do nothing
  } else {
    if(!allowNA && any(is.na(arg)))
      stop("Argument '", name, "' must not be NA or NaN.", call.=FALSE)
    if(!is.numeric(arg))
      stop("Argument '", name, "' must be numeric.", call.=FALSE)
    if(any(arg < min | arg > max, na.rm=TRUE))
      stop("Argument '", name, "' must be between ", min, " and ", max, ".", call.=FALSE)
  }
}

stopifNegative <- function(arg, allowNA=FALSE, allowZero=TRUE) {
  name <- deparse(substitute(arg))
  if(allowNA && all(is.na(arg))) {  # An all-NA vector is logical, but ok.
    # do nothing
  } else {
    if(!allowNA && any(is.na(arg)))
      stop("Argument '", name, "' must not be NA or NaN.", call.=FALSE)
    if(!is.numeric(arg))
      stop("Argument '", name, "' must be numeric.", call.=FALSE)
    if(allowZero) {
      if(any(arg < 0, na.rm=TRUE))
        stop("Argument '", name, "' must be non-negative.", call.=FALSE)
    } else {
      if(any(arg <= 0, na.rm=TRUE))
        stop("Argument '", name, "' must be greater than 0.", call.=FALSE)
    }
  }
}

