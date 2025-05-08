#' Predicted values for dose-response network meta-analysis
#'
#' @description
#' This function provides the predicted values based on the results of dose-response
#' network meta-analysis.
#'
#' @param object An object of class netdose (mandatory).
#' @param agent1 An optional character string specifying the first agent
#'   to be used for the prediction. By default, all agents are used.
#' @param dose1 An optional numeric vector specifying custom doses for the
#'   prediction. By default, the doses are set to the common observed doses as defined in the data.
#' @param agent2 An optional character string specifying the second agent
#'   to be used for the prediction. By default, the reference agent is used.
#' @param dose2 An optional numeric vector specifying the dose for the second
#'   agent. By default the common dose of the second agent is used.
#' @param ... Additional arguments (ignored).
#'
#' @details
#' The predict.netdose function calculates predicted effects for specified doses of one or more
#' agents, based on a dose-response network meta-analysis.
#' It supports both linear and non-linear dose-response relationships, accommodating various
#' modeling methods including linear, exponential, fractional polynomials, restricted cubic splines (RCS),
#' and quadratic relationships.
#' This function is particularly useful for exploring comparative effectiveness
#' at specific dose levels of the agents, facilitating the interpretation of
#' complex dose-response relationships in a network meta-analysis setting.
#' By allowing predictions for multiple combinations of agents and doses,
#' it offers flexibility in evaluating hypothetical scenarios or estimating effects
#' for doses outside the directly observed range (where extrapolation is appropriate).
#'
#' @return
#' A data frame with additional class \code{predict.netdose} containing the
#' following variables:
#' \item{agent1, dose1, agent2, dose2}{As defined above}
#' \item{pred}{A numeric vector with the predicted effects}
#' \item{se.pred}{A numeric vector with standard errors of the predicted effects}
#' \item{lower}{A numeric vector specifying the lower bounds of the predicted values}
#' \item{upper}{A numeric vector specifying the upper bounds of the predicted values}
#'

#' @author Maria Petropoulou <maria.petropoulou@@uniklinik-freiburg.de>,
#'   Guido Schwarzer <guido.schwarzer@@uniklinik-freiburg.de>
#'
#' @references
#' Petropoulou et al. (2025):
#' Network meta-analysis with dose-response relationships.
#'
#' @examples
#' # Use a subset of 5 studies from anesthesia data
#' anesthesia_subset <- subset(anesthesia, study %in% unique(anesthesia$study)[1:5])
#'
#' # Prepare data for DR-NMA
#' dat <- pairwise(
#'   agent = list(agent1, agent2, agent3),
#'   event = list(event1, event2, event3),
#'   n = list(n1, n2, n3),
#'   dose = list(dose1, dose2, dose3),
#'   data = anesthesia_subset,
#'   studlab = study,
#'   append = FALSE
#' )
#'
#'
#' # Perform DR-NMA with a linear dose-response function
#' dr1 <- netdose(
#'   TE, seTE, agent1, dose1, agent2,
#'   dose2, studlab,
#'   data = dat
#' )
#'
#' # Predicted values
#' pred1 <- predict(dr1)
#' 
#' @method predict netdose
#' @export


predict.netdose <- function(object,
                            agent1 = NULL, dose1 = NULL,
                            agent2 = NULL, dose2 = NULL,
                            ...) {
  #
  # Check class
  #
  chkclass(object, "netdose")

  #
  #  Check arguments
  #
  if (!is.null(agent1)) {
    agent1 <- setchar(agent1, object$agents)
  } else {
    agent1 <- object$agents[!(object$agents %in% object$inactive)]
  }
  #
  if (!is.null(agent2)) {
    agent2 <- setchar(agent2, object$agents)
  } else {
    agent2 <- object$reference.group
  }
  #
  if (!is.null(dose1)) {
    chknumeric(dose1)
  } else {
    dose1 <- object$common.dose[agent1]
  }
  #
  if (!is.null(dose2)) {
    chknumeric(dose2)
  } else {
    dose2 <- object$common.dose[agent2]
  }
  #
  k.Pred <- max(c(length(agent1), length(dose1), length(agent2), length(dose2)))
  #
  if (k.Pred > 1) {
    if (length(agent1) == 1) {
      agent1 <- rep(agent1, k.Pred)
    } else {
      chklength(agent1, k.Pred,
        text =
          paste0(
            "Input for argument 'agent1' must be of length ",
            k.Pred, "."
          )
      )
    }
    #
    if (length(dose1) == 1) {
      dose1 <- rep(dose1, k.Pred)
    } else {
      chklength(dose1, k.Pred,
        text =
          paste0(
            "Input for argument 'dose1' must be of length ",
            k.Pred, "."
          )
      )
    }
    #
    if (length(agent2) == 1) {
      agent2 <- rep(agent2, k.Pred)
    } else {
      chklength(agent2, k.Pred,
        text =
          paste0(
            "Input for argument 'agent2' must be of length ",
            k.Pred, "."
          )
      )
    }
    #
    if (length(dose2) == 1) {
      dose2 <- rep(dose2, k.Pred)
    } else {
      chklength(dose2, k.Pred,
        text =
          paste0(
            "Input for argument 'dose2' must be of length ",
            k.Pred, "."
          )
      )
    }
  }
  #
  if (object$random) {
    # Coefficients
    coef <- object$TE.drnma.random
    se.coef <- object$seTE.drnma.random
    lower.coef <- object$lower.drnma.random
    upper.coef <- object$upper.drnma.random
  } else {
    # Coefficients
    coef <- object$TE.drnma.common
    se.coef <- object$seTE.drnma.common
    lower.coef <- object$lower.drnma.common
    upper.coef <- object$upper.drnma.common
  }

  #
  # Loop through requested predictions
  #
  pred <- se.pred <- lower <- upper <- vector("numeric", k.Pred)
  #
  for (i in seq_len(k.Pred)) {
    if (object$method %in% c("linear", "exponential", "fp1")) {
      if (object$method == "linear") {
        g1 <- dose2dose
        #
        param1 <- NULL
      } else if (object$method == "exponential") {
        g1 <- dose2exp
        #
        param1 <- NULL
      } else if (object$method == "fp1") {
        g1 <- dose2fp
        #
        param1 <- object$param
      }
      #
      pred[i] <-
        g1(dose1[i], param1) * sel_coef(coef, agent1[i]) -
        g1(dose2[i], param1) * sel_coef(coef, agent2[i])
      #
      se.pred[i] <-
        sqrt((g1(dose1[i], param1) * sel_coef(se.coef, agent1[i]))^2 +
          (g1(dose2[i], param1) * sel_coef(se.coef, agent2[i]))^2)
      #
      lower[i] <-
        g1(dose1[i], param1) * sel_coef(lower.coef, agent1[i]) -
        g1(dose2[i], param1) * sel_coef(lower.coef, agent2[i])
      #
      upper[i] <-
        g1(dose1[i], param1) * sel_coef(upper.coef, agent1[i]) -
        g1(dose2[i], param1) * sel_coef(upper.coef, agent2[i])
    } else {
      if (object$method == "quadratic") {
        g1 <- dose2dose
        g2 <- dose2poly
        #
        param1 <- NULL
        param2 <- NULL
      } else if (object$method == "fp2") {
        g1 <- dose2fp
        g2 <- dose2fp
        #
        param1 <- object$param[1]
        param2 <- object$param[2]
      } else if (object$method == "rcs") {
        g1 <- dose2dose
        g2 <- dose2rcs
        #
        param1 <- NULL
        param2 <- object$param
      }
      #
      pred[i] <-
        g1(dose1[i], param1) * sel_coef(coef, agent1[i]) +
        g2(dose1[i], param2) * sel_coef(coef, agent1[i], 2) -
        g1(dose2[i], param1) * sel_coef(coef, agent2[i]) +
        g2(dose2[i], param2) * sel_coef(coef, agent2[i], 2)
      #
      se.pred[i] <-
        sqrt((g1(dose1[i], param1) * sel_coef(se.coef, agent1[i]))^2 +
          (g2(dose1[i], param2) * sel_coef(se.coef, agent1[i], 2))^2 +
          (g1(dose2[i], param1) * sel_coef(se.coef, agent2[i]))^2 +
          (g2(dose2[i], param2) * sel_coef(se.coef, agent2[i], 2))^2)
      #
      lower[i] <-
        g1(dose1[i], param1) * sel_coef(lower.coef, agent1[i]) +
        g2(dose1[i], param2) * sel_coef(lower.coef, agent1[i], 2) -
        g1(dose2[i], param1) * sel_coef(lower.coef, agent2[i]) +
        g2(dose2[i], param2) * sel_coef(lower.coef, agent2[i], 2)
      #
      upper[i] <-
        g1(dose1[i], param1) * sel_coef(upper.coef, agent1[i]) +
        g2(dose1[i], param2) * sel_coef(upper.coef, agent1[i], 2) -
        g1(dose2[i], param1) * sel_coef(upper.coef, agent2[i]) +
        g2(dose2[i], param2) * sel_coef(upper.coef, agent2[i], 2)
    }
  }


  res <- data.frame(agent1, dose1, agent2, dose2, pred, se.pred, lower, upper)
  #
  attr(res, "level") <- object$level
  attr(res, "sm") <- object$sm
  attr(res, "backtransf") <- object$backtransf
  #
  class(res) <- c("predict.netdose", class(res))
  #
  res
}
