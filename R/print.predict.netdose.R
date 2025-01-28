#' Print method for objects of class predict
#'
#' @param x An object of class \code{predict.netdose}.
#' @param backtransf A logical indicating whether printed results
#'   should be back transformed. If \code{backtransf=TRUE}, results
#'   for \code{sm="OR"} are printed as odds ratios rather than log
#'   odds ratios, for example.
#' @param digits Minimal number of significant digits, see
#'   \code{print.default}.
#' @param big.mark A character used as thousands separator.
#' @param \dots Additional arguments (ignored).
#'
#' @method print predict.netdose
#' @export
#'
#'


print.predict.netdose <- function(x,
                                  backtransf = attr(x, "backtransf"),
                                  digits = gs("digits"),
                                  big.mark = gs("big.mark"),
                                  ...) {
  # Check class
  chkclass(x, "predict.netdose")
  
  # Check arguments
  chklogical(backtransf)
  chknumeric(digits)
  #
  sm <- attr(x, "sm")
  smlab <- sm
  #
  cilab <- paste0(round(100 * attr(x, "level"), 1), "%-CI")
  
  #
  # Back-transform results
  #
  if (is_relative_effect(sm) & backtransf) {
    x$pred <- exp(x$pred)
    x$lower <- exp(x$lower)
    x$upper <- exp(x$upper)
  }
  else if (is_relative_effect(sm) & !backtransf)
    smlab <- paste0("log", sm)
  
  #
  # Round results
  #
  x$pred <- round(x$pred, digits = digits)
  x$lower <- round(x$lower, digits = digits)
  x$upper <- round(x$upper, digits = digits)
  
  #
  # Print predictions
  #
  res <- x
  #
  res$comparison <-
    paste(paste0("'", res$agent1, " ", res$dose1, "'"),
          " vs ",
          paste0("'", res$agent2, " ", res$dose2, "'"))
  #
  res$ci <-
    formatCI(
      formatN(x$lower, digits = digits, big.mark = big.mark),
      formatN(x$upper, digits = digits, big.mark = big.mark))
  #
  res <- res[, c("comparison", "pred", "ci")]
  names(res) <- c("comparison", smlab, cilab)
  #
  names(res)[names(res) == "ci"] <-
    
  #
  prmatrix(res, rowlab = rep("", nrow(res)), quote = FALSE, right = TRUE)
  
  invisible(NULL)
}
