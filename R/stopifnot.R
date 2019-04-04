
# Sanity checks

# Not exported.

# These functions are very basic but simplify input checks while giving
#   informative error messages.
# Without checks, users may get mysterious error messages, eg,
#   "is.na() applied to non-(list or vector) of type 'closure'".
# base::stopifnot() is better, but error messages are still abstruse, eg,
#  "p >= 0 & p <= 1 is not TRUE".

stopifnotNumeric <- function(arg, allowNA=FALSE) {
  name <- deparse(substitute(arg))
  if(allowNA && all(is.na(arg))) {
    # do nothing
  } else {
    if(!allowNA && any(is.na(arg)))
      stop("Argument '", name, "' must not contain NA or NaN.", call.=FALSE)
    if(!is.numeric(arg))
      stop("Argument '", name, "' must be numeric.", call.=FALSE)
  }
}

stopifnotEqualLength <- function(arg1, arg2) {
  name1 <- deparse(substitute(arg1))
  name2 <- deparse(substitute(arg2))
  if(length(arg1) != length(arg2))
   stop("Lengths of arguments '", name1, "' and '", name2, "' must be equal.", call.=FALSE)
}

stopifnotGreaterthan <- function(arg, value, allowNA=FALSE) {
  name <- deparse(substitute(arg))
  if(allowNA && all(is.na(arg))) {
    # do nothing
  } else {
    if(!allowNA && any(is.na(arg)))
      stop("Argument '", name, "' must not contain NA or NaN.", call.=FALSE)
    if(!is.numeric(arg))
      stop("Argument '", name, "' must be numeric.", call.=FALSE)
    if(any(arg <= value))
      if(allowNA) {
        stop("Argument '", name, "' must be greater than ", value, ", or NA.", call.=FALSE)
      } else {
        stop("Argument '", name, "' must be greater than ", value, ".", call.=FALSE)
      }
  }
}

stopifnotLessthan <- function(arg, value, allowNA=FALSE) {
  name <- deparse(substitute(arg))
  if(allowNA && all(is.na(arg))) {
    # do nothing
  } else {
    if(!allowNA && any(is.na(arg)))
      stop("Argument '", name, "' must not contain NA or NaN.", call.=FALSE)
    if(!is.numeric(arg))
      stop("Argument '", name, "' must be numeric.", call.=FALSE)
    if(any(arg >= value)) {
      if(allowNA) {
        stop("Argument '", name, "' must be less than ", value, ", or NA.", call.=FALSE)
      } else {
        stop("Argument '", name, "' must be less than ", value, ".", call.=FALSE)
      }
    }
  }
}

stopifnotInteger <- function(arg, allowNA=FALSE) {
  name <- deparse(substitute(arg))
  if(allowNA && all(is.na(arg))) {
    # do nothing
  } else {
    if(!allowNA && any(is.na(arg)))
      stop("Argument '", name, "' must not contain NA or NaN.", call.=FALSE)
    if(!is.numeric(arg))
      stop("Argument '", name, "' must be numeric.", call.=FALSE)
    if(!all(arg%%1 == 0)) {
      if(allowNA) {
        stop("Argument '", name, "' must be integer (whole number), or NA.", call.=FALSE)
      } else {
        stop("Argument '", name, "' must be integer (whole number).", call.=FALSE)
      }
    }
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

stopifnotLength <- function(arg, length, allow1=FALSE) {
  name <- deparse(substitute(arg))
  if(allow1 && length(arg) == 1) {
    # do nothing
  } else {
    if(length(arg) != length) {
      if(allow1) {
        stop("Argument '", name, "' must have length 1 or ", length, ".", call.=FALSE)
      } else {
        stop("Argument '", name, "' must have length ", length, ".", call.=FALSE)
      }
    }
  }
}

stopifnotProbability <- function(arg, allowNA=FALSE) {
  name <- deparse(substitute(arg))
  if(allowNA && all(is.na(arg))) {  # An all-NA vector is logical, but ok.
    # do nothing
  } else {
    if(!allowNA && any(is.na(arg)))
      stop("Argument '", name, "' must not contain NA or NaN.", call.=FALSE)
    if(!is.numeric(arg))
      stop("Argument '", name, "' must be numeric.", call.=FALSE)
    if(any(arg < 0 | arg > 1, na.rm=TRUE)) {
      if(allowNA) {
        stop("Argument '", name, "' must be a probability between 0 and 1, or NA.", call.=FALSE)
      } else {
        stop("Argument '", name, "' must be a probability between 0 and 1.", call.=FALSE)
      }
    }
  }
}

stopifnotBetween <- function(arg, min, max, allowNA=FALSE) {
  name <- deparse(substitute(arg))
  if(allowNA && all(is.na(arg))) {  # An all-NA vector is logical, but ok.
    # do nothing
  } else {
    if(!allowNA && any(is.na(arg)))
      stop("Argument '", name, "' must not contain NA or NaN.", call.=FALSE)
    if(!is.numeric(arg))
      stop("Argument '", name, "' must be numeric.", call.=FALSE)
    if(any(arg < min | arg > max, na.rm=TRUE)) {
      if(allowNA) {
        stop("Argument '", name, "' must be between ", min, " and ", max, ", or NA.", call.=FALSE)
      } else{
        stop("Argument '", name, "' must be between ", min, " and ", max, ".", call.=FALSE)
      }
    }
  }
}

stopifNegative <- function(arg, allowNA=FALSE, allowZero=TRUE) {
  name <- deparse(substitute(arg))
  if(allowNA && all(is.na(arg))) {  # An all-NA vector is logical, but ok.
    # do nothing
  } else {
    if(!allowNA && any(is.na(arg)))
      stop("Argument '", name, "' must not contain NA or NaN.", call.=FALSE)
    if(!is.numeric(arg))
      stop("Argument '", name, "' must be numeric.", call.=FALSE)
    if(allowZero) {
      if(any(arg < 0, na.rm=TRUE)) {
        if(allowNA) {
          stop("Argument '", name, "' must be non-negative, or NA.", call.=FALSE)
        } else {
          stop("Argument '", name, "' must be non-negative.", call.=FALSE)
        }
      }
    } else {
      if(any(arg <= 0, na.rm=TRUE)) {
        if(allowNA) {
          stop("Argument '", name, "' must be greater than 0, or NA.", call.=FALSE)
        } else {
          stop("Argument '", name, "' must be greater than 0.", call.=FALSE)
        }
      }
    }
  }
}

