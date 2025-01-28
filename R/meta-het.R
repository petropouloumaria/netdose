#
# Auxiliary functions to calculate heterogeneity measures (from R package meta)
#
# Package: netdose
# Author: Guido Schwarzer <guido.schwarzer@@uniklinik-freiburg.de>
# License: GPL (>= 2)
#

calcH <- function(Q, df, level) {
  #
  # Calculate H
  # Higgins & Thompson (2002), Statistics in Medicine, 21, 1539-58
  #
  k <- df + 1
  #
  if (!is.na(k)) {
    if (k > 1)
      H <- sqrt(Q / (k - 1))
    else
      H <- NA
    #
    selogH <- ifelse(Q > k,
              ifelse(k >= 2,
                     0.5 * (log(Q) - log(k - 1)) /
                     (sqrt(2 * Q) - sqrt(2 * k - 3)),
                     NA),
              ifelse(k > 2,
                     sqrt(1 / (2 * (k - 2)) * (1 - 1 / (3 * (k - 2)^2))),
                     NA))
  }
  else {
    H <- NA
    selogH <- NA
  }
  #
  tres <- ci(log(max(c(H, 1))), selogH, level)
  #
  res <- list(TE = exp(tres$TE),
              lower = max(exp(tres$lower), 1),
              upper = max(exp(tres$upper), 1))
  
  res
}

isquared <- function(Q, df, level) {
  #
  # Calculate I-Squared
  # Higgins & Thompson (2002), Statistics in Medicine, 21, 1539-58
  #
  tres <- calcH(Q, df, level)
  #
  func.t <- function(x) (x^2 - 1) / x^2
  #
  res <- list(TE = func.t(tres$TE),
              lower = func.t(tres$lower),
              upper = func.t(tres$upper))
  
  res
}

pvalQ <- function(Q, df, lower.tail = FALSE) {
  #
  if (length(df) == 1 & length(Q) > 1)
    df <- rep(df, length(Q))
  else if (length(df) > 1 & length(Q) == 1)
    Q <- rep(Q, length(df))
  else if (length(df) > 1 & length(Q) > 1 & length(df) != length(Q))
    stop("Length of arguments 'Q' and 'df' do not match.")
  #
  res <- ifelse(is.na(df) | df < 1,
                NA,
                pchisq(Q, df, lower.tail = lower.tail))
  #
  res
}
