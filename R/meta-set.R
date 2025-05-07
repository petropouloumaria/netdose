#
# Auxiliary functions to set arguments (from R package meta)
#
# Package: netdose
# Author: Guido Schwarzer <guido.schwarzer@uniklinik-freiburg.de>
# License: GPL (>= 2)
#

setchar <- function(x, val, text, list = FALSE, name = NULL,
                    stop.at.error = TRUE, addtext = "",
                    return.NULL = TRUE, nchar.equal = FALSE,
                    setNA = FALSE) {
  val <- unique(val)
  #
  if (is.null(name))
    name <- deparse(substitute(x))
  nval <- length(val)
  #
  if (is.numeric(x)) {
    numeric.x <- TRUE
    idx <- x
    idx[idx < 1] <- NA
    idx[idx >= nval + 1] <- NA
  }
  else {
    numeric.x <- FALSE
    #
    if (length(unique(tolower(x))) != length(unique(x)) |
        length(unique(tolower(val))) != length(unique(val)))
      idx <- charmatch(x, val, nomatch = NA)
    else
      idx <- charmatch(tolower(x), tolower(val), nomatch = NA)
  }
  #
  if ((anyNA(idx) || any(idx == 0)) && !setNA) {
    if (list)
      first <- "List element '"
    else
      first <- "Argument '"
    #
    if (missing(text)) {
      if (numeric.x) {
        if (nval == 1)
          vlist <- "1"
        else if (nval == 2)
          vlist <- "1 or 2"
        else
          vlist <- paste("between 1 and", nval)
      }
      else {
        if (nval == 1)
          vlist <- paste0('"', val, '"')
        else if (nval == 2)
          vlist <- paste0('"', val, '"', collapse = " or ")
        else
          vlist <- paste0(paste0('"', val[-nval], '"', collapse = ", "),
                          ', or ', '"', val[nval], '"')
      }
      #
      if (stop.at.error)
        stop(first, name, "' must be ", vlist, addtext, ".", call. = FALSE)
      else {
        if (return.NULL)
          return(NULL)
        else
          return(x)
      }
    }
    else {
      if (stop.at.error)
        stop(first, name, "' ", text, ".", call. = FALSE)
      else {
        if (return.NULL)
          return(NULL)
        else
          return(x)
      }
    }
  }
  #
  if (is.null(x))
    return(NULL)
  else
    res <- val[idx]
  #
  if (nchar.equal && nchar(res) != nchar(x))
    res <- x
  #
  res
}

setstudlab <- function(x, k) {
  #
  # Set study labels
  #
  if (is.null(x)) {
    if (k == 1)
      x <- ""
    else
      x <- seq(along = rep(1, k))
  }
  #
  if (is.factor(x))
    x <- as.character(x)
  #
  x
}

setVal <- function(data, varname, default = NULL) {
  if (isCol(data, varname))
    return(data[[varname]])
  else
    return(default)
}

setsort <- function(sort, n, text) {
  if (is.null(sort))
    res <- seq_len(n)
  else {
    chklength(sort, n,
              text = paste0("Argument '", deparse(substitute(sort)),
                           "' must be of same length as ",
                           "number of ", text, "."))
    #
    res <- sort
    if (!(is.numeric(res) & min(res) == 1 & max(res) == n))
      res <- order(res)
  }
  #
  res
}
