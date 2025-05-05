#'  Dose-response curve plot
#'
#' @description
#' Generates a dose-response plot based on the results of a dose-response
#' network meta-analysis (DR-NMA). The plot visualizes predicted dose-response
#' curves alongside observed responses for easy interpretation of model outputs.
#'
#' @param x An object of class netdose (mandatory).
#' @param pooled A character string indicating whether results for the
#'   common (\code{"common"}) or random effects model
#'   (\code{"random"}) should be plotted.  Abbreviations are allowed.
#'   Defaults to "random" if the input object specifies a random effects model;
#'   otherwise, defaults to "common".
#' @param only.direct A logical value indicating whether only the study results
#'   of direct comparisons with the reference agent for the observed data should
#'   be shown in the plot. Defaults to \code{TRUE}.
#' @param col.direct The color used for points representing direct comparisons.
#'   By default, \code{"black"} when \code{only.direct = TRUE}; otherwise, \code{"green"}.
#' @param col.indirect The color used for points representing indirect comparisons.
#'   Defaults to \code{"red"}.
#' @param \dots Additional arguments. Currently ignored, but included for potential
#'   future extensions or compatibility with generic plotting functions.
#'
#'
#' @details
#' The function plots the dose-response curve alongside the observed responses:
#' - The vertical axis represents the dose range, which is defined from 0 to the
#'   maximum observed dose, with 100 evenly spaced points generated within this range.
#' - The horizontal axis represents the predicted response values, calculated using
#'   the \code{predict.netdose} function.
#'
#' The plot includes shaded confidence intervals for the predicted dose-response
#' curve. Observed responses are overlaid for comparison, differentiated into
#' direct and indirect comparisons with customizable colors.
#'
#' This visualization aids in understanding the relationship between dose and
#' response while validating the model fit against observed data.
#'
#' @keywords hplot
#'
#' @author Maria Petropoulou <maria.petropoulou@uniklinik-freiburg.de>,
#'  Guido Schwarzer <guido.schwarzer@uniklinik-freiburg.de>
#'
#' @examples
#' # Use a subset of 8 studies
#' anesthesia_subset <-
#'   subset(anesthesia, study %in% unique(anesthesia$study)[1:20])
#'
#' # Prepare data for DR-NMA
#' dat <- pairwise(
#'   agent = list(agent1, agent2, agent3, agent4, agent5),
#'   event = list(event1, event2, event3, event4, event5),
#'   n = list(n1, n2, n3, n4, n5),
#'   dose = list(dose1, dose2, dose3, dose4, dose5),
#'   data = anesthesia_subset,
#'   studlab = study,
#'   append = FALSE
#' )
#'
#' # DR-NMA with linear dose-response function
#' dr1 <- netdose(TE, seTE, agent1, dose1, agent2,
#'   dose2, studlab,
#'   data = dat
#' )
#'
#' # Dose-response plot
#' plot(dr1)
#'
#' @method plot netdose
#' @export
#'


plot.netdose <- function(x, pooled = if (x$random) "random" else "common",
                         only.direct = TRUE,
                         col.direct = if (only.direct) "black" else "green",
                         col.indirect = "red",
                         ...) {
  # Check class
  chkclass(x, "netdose")
  #
  data <- x$data
  #
  pooled <- setchar(pooled, c("common", "random"))
  chklogical(only.direct)

  # Get rid of warnings "no visible binding for global variable"
  .agent1 <- .agent2 <- dose1 <- TE.adj <-
    pred <- lower.pred <- upper.pred <- TE <- NULL

  # Coefficients
  if (pooled == "random") {
    coef <- x$TE.drnma.random
  } else {
    coef <- x$TE.drnma.common
  }

  # exclude reference group from the plot
  coef <- coef[!names(coef) %in% x$reference.group]

  # Create a list to store plots
  plots_list <- list()

  # Loop through each agent and create a plot
  for (i in names(coef)) {
    dat.i <- subset(data, .agent1 == i | .agent2 == i)
    # Drop comparisons using the same agent
    dat.i <- subset(dat.i, !(.agent1 == i & .agent2 == i))
    #
    wo <- dat.i$.agent1 != i & dat.i$.agent2 == i
    #
    if (any(wo)) {
      dat.i$.TE[wo] <- -dat.i$.TE[wo]
      #
      tagent1 <- dat.i$.agent1
      dat.i$.agent1[wo] <- dat.i$.agent2[wo]
      dat.i$.agent2[wo] <- tagent1[wo]
      #
      tdose1 <- dat.i$.dose1
      dat.i$.dose1[wo] <- dat.i$.dose2[wo]
      dat.i$.dose2[wo] <- tdose1[wo]
    }
    # Active vs active
    active <- !(dat.i$.agent2 %in% x$inactive | dat.i$dose2 == 1)
    #
    pred2.i <- predict(x, agent1 = dat.i$.agent2, dose1 = dat.i$.dose2)
    #
    # Data set with observed intervention effects
    dat.obs <- data.frame(
      agent1 = dat.i$.agent1, dose1 = dat.i$.dose1,
      agent2 = dat.i$.agent2, dose2 = dat.i$.dose2,
      TE = dat.i$.TE
    )
    # Add predicted effect of second intervention
    if (!only.direct) {
      dat.obs$TE.adj <- ifelse(active, dat.obs$TE + pred2.i$pred, dat.obs$TE)
    } else {
      dat.obs$TE.adj <- dat.obs$TE
      dat.obs <- dat.obs[!active, , drop = FALSE]
    }

    # Get the line with predicted values
    seq1 <- seq(0, max(dat.i$.dose1, na.rm = TRUE), length.out = 101)[-1]
    predline.i <- predict(x, agent1 = i, dose1 = seq1)
    #
    linedata <- data.frame(
      dose1 = predline.i$dose1,
      pred = predline.i$pred,
      lower.pred = predline.i$lower,
      upper.pred = predline.i$upper
    )

    # Create plot for the current agent
    plot <- ggplot(linedata, aes(x = dose1, y = pred)) +
      geom_ribbon(
        aes(
          ymin = lower.pred,
          ymax = upper.pred,
        ),
        fill = "darkgrey", alpha = 0.2
      ) +
      geom_line(color = "blue") +
      labs(y = NULL, x = NULL, title = i) +
      theme(plot.title = element_text(size = 10)) +
      theme(
        panel.background = element_rect(fill = "transparent", color = NA), # Set panel background to transparent
        plot.background = element_rect(fill = "transparent", color = NA)
      )
    #
    # Add observed values
    #
    if (only.direct) {
      plot <- plot +
        geom_point(
          data = dat.obs, aes(y = TE.adj), color = col.direct,
          size = 1
        )
    } else {
      plot <- plot +
        geom_point(
          data = dat.obs, aes(y = TE.adj),
          color = ifelse(active, col.indirect, col.direct), size = 1
        )
    }

    plots_list[[i]] <- plot
  } # end for loop

  #
  # Arrange multiple plots in a single figure
  #
  if (x$method == "fp1") {
    if (is.null(x$param)) {
      x$param <- -0.5
    }
    x$method <- paste0("FP1 (p=", x$param, ")")
  } else if (x$method == "fp2") {
    if (is.null(x$param)) {
      x$param <- c(-0.5, 1)
    }
    #
    x$method <- paste0("FP2 (p1=", x$param[1], ", p2=", x$param[2], ")")
  } else if (x$method == "rcs") {
    x$method <- paste0("RCS (p=", paste(x$param, collapse = ", "), ")")
  } else {
    x$method <-
      paste0(toupper(substring(x$method, 1, 1)), substring(x$method, 2))
  }
  #
  grid.arrange(
    grobs = plots_list, ncol = 7,
    top = paste(x$method, "dose-response plot"),
    left = xlab(x$sm, FALSE),
    bottom = "Dose"
  )
}
