setref <- function(reference.group, levs, length = 1,
                   varname = "reference.group", error.text) {
  
  if (missing(error.text)) {
    text.start <- paste0("Argument '", varname, "'")
    text.within <- paste0("argument '", varname, "'")
  }
  else {
    text.start <- paste0(toupper(substring(error.text, 1, 1)),
                         substring(error.text, 2))
    text.within <- error.text
  }
  
  
  if (length && length(reference.group) != length)
    stop(text.start,
         if (length == 1)
           " must be a numeric or a character string"
         else
           paste(" must be a numeric of character vector of length", length),
         ".",
         call. = FALSE)
  ##
  if (is.numeric(reference.group)) {
    if (any(is.na(reference.group)))
      stop("Missing value not allowed in ", text.within, ".",
           call. = FALSE)
    if (!all(reference.group %in% seq_len(length(levs))))
      stop(paste0(text.start, " must ",
                  if (length == 1) "be any of the " else "contain ",
                  "integers from 1 to ",
                  length(levs), "."),
           call. = FALSE)
    res <- levs[reference.group]
  }
  else if (is.character(reference.group)) {
    if (any(is.na(reference.group)))
      stop("Missing value not allowed in ", text.within, ".",
           call. = FALSE)
    ##
    if (length(unique(levs)) == length(unique(tolower(levs))))
      idx <- charmatch(tolower(reference.group), tolower(levs), nomatch = NA)
    else {
      idx1 <- charmatch(reference.group, levs, nomatch = NA)
      idx2 <- charmatch(tolower(reference.group), tolower(levs), nomatch = NA)
      if (anyNA(idx1) & !anyNA(idx2))
        idx <- idx2
      else
        idx <- idx1
    }
    ##
    if (anyNA(idx) || any(idx == 0))
      stop("Admissible values for ", text.within, ":\n  ",
           paste(paste0("'", levs, "'"), collapse = " - "),
           "\n  (unmatched value", if (sum(is.na(idx)) > 1) "s",
           ": ",
           paste(paste0("'", reference.group[is.na(idx)], "'"),
                 collapse = " - "),
           ")",
           call. = FALSE)
    res <- levs[idx]
  }
  
  res
}
