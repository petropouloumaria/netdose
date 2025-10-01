#' Dose-response curve plot
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
#'   By default, \code{"black"} when \code{only.direct = TRUE}; otherwise,
#'   \code{"green"}.
#' @param col.indirect The color used for points representing indirect
#'   comparisons. Defaults to \code{"red"}.
#' @param agents Optional character vector specifying which agents to include
#'   in the plot. If NULL, all agents will be plotted.
#' @param ylim Optional numeric vector of length 2 specifying the y-axis limits.
#'   If NULL, limits are determined automatically.
#' @param benchmark.threshold Numeric; benchmark response level (e.g., 0.1 for
#'   10 percent). Used to compute Benchmark Dose Lower Confidence Limit (BMDL).
#' @param plateau.threshold Numeric; threshold for identifying the plateau in
#'   the dose-response curve. Defines the minimum absolute change in predicted
#'   response between adjacent dose levels, below which the response is
#'   considered stable (i.e., plateau has been reached). Used to calculate the
#'   Maximum Effective Dose (MED). Default: \code{0.0001}.
#' @param col.line Colour for the dose-response line.
#' @param col.bmdl Colour for the BMLD line.
#' @param col.med Colour for the MED line.
#' @param legend A logical value indicating whether to print a legend.
#' @param \dots Additional arguments. Currently ignored, but included for
#'   potential future extensions or compatibility with generic plotting
#'   functions.
#'
#' @details
#' The function plots the dose-response curve alongside the observed responses:
#' 
#' \itemize{
#' \item The vertical axis represents the dose range, which is defined from 0
#'   to the maximum observed dose, with 100 evenly spaced points generated
#'   within this range.
#' \item The horizontal axis represents the predicted response values,
#'   calculated using the \code{predict.netdose} function.
#' }
#' 
#' The plot includes shaded confidence intervals for the predicted dose-response
#' curve. Observed responses are overlaid for comparison, differentiated into
#' direct and indirect comparisons with customizable colors.
#' 
#' The function also optionally displays the Benchmark Dose (BMD) and the
#' Benchmark Dose Lower Confidence Limit (BMDL), based on a user-defined
#' benchmark response threshold (e.g., 0.01 for 10 percent increase).
#'
#' If the model indicates that the predicted response stabilizes beyond a
#' certain dose level, the function estimates and plots the Maximum Effective
#' Dose (MED) — the smallest dose beyond which the predicted response increases
#' less than a given threshold (controlled via \code{plateau.threshold}). MED
#' is shown only if it occurs after the BMDL, ensuring biological and
#' statistical coherence.
#' 
#' @keywords hplot
#'
#' @author Maria Petropoulou <maria.petropoulou@@uniklinik-freiburg.de>,
#'  Guido Schwarzer <guido.schwarzer@@uniklinik-freiburg.de>
#'
#' @examples
#' # Use a subset of 5 studies from anesthesia data
#' anesthesia_subset <-
#'   subset(anesthesia, study %in% unique(anesthesia$study)[1:5])
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
#' # DR-NMA with linear dose-response function
#' dr1 <- netdose(TE, seTE, agent1, dose1, agent2,
#'   dose2, studlab,
#'   data = dat
#' )
#'
#' # Dose-response plot
#' plot(dr1)
#' 
#' @return No return value, called for side effects (generates a plot).
#' @method plot netdose
#' @export

plot.netdose <- function(x, pooled = if (x$random) "random" else "common",
                         only.direct = TRUE,
                         col.direct = if (only.direct) "black" else "green",
                         col.indirect = "red",
                         agents = NULL,
                         ylim = NULL,
                         benchmark.threshold = NULL,
                         plateau.threshold = NULL,
                         col.line = "blue",
                         col.bmdl = "purple",
                         col.med = "gray40",
                         legend = !only.direct,
                         ...) {
  # Check class
  chkclass(x, "netdose")
  #
  data <- x$data
  #
  pooled <- setchar(pooled, c("common", "random"))
  chklogical(only.direct)
  #
  if (!is.null(benchmark.threshold))
    chknumeric(benchmark.threshold, min = 0, zero = TRUE, length = 1)
  #
  if (!is.null(plateau.threshold))
    chknumeric(plateau.threshold, min = 0, zero = TRUE, length = 1)
  #
  if (length(col.line) != 1)
    stop("Argument 'col.line' must be a single color.", call. = FALSE)
  #
  if (length(col.bmdl) != 1)
    stop("Argument 'col.bmdl' must be a single color.", call. = FALSE)
  #
  if (length(col.med) != 1)
    stop("Argument 'col.med' must be a single color.", call. = FALSE)
  #
  chklogical(legend)
  
  # Get rid of warnings "no visible binding for global variable"
  .agent1 <- .agent2 <- dose1 <- TE.adj <-
    pred <- lower.pred <- upper.pred <- TE <-
    bmdl <- med <- type <- NULL
  
  # Coefficients
  if (pooled == "random") {
    coef <- x$TE.drnma.random
  }
  else {
    coef <- x$TE.drnma.common
  }
  
  # exclude reference group from the plot
  coef <- coef[!names(coef) %in% x$reference.group]
  
  if (!is.null(agents)) {
    coef <- coef[names(coef) %in% agents]
  }
  
  
  dat <- preddat <- bmdldat <- meddat <- data.frame()
  #
  # Loop through each agent and create data sets
  #
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
    active <- !(dat.i$.agent2 %in% x$inactive | dat.i$.dose2 == 0)
    #
    pred2.i <- predict(x, agent1 = dat.i$.agent2, dose1 = dat.i$.dose2)
    #
    # Data set with observed intervention effects
    obsdat.i <- data.frame(
      agent = i,
      agent1 = dat.i$.agent1, dose1 = dat.i$.dose1,
      agent2 = dat.i$.agent2, dose2 = dat.i$.dose2,
      TE = dat.i$.TE
    ) %>%
      mutate(type = ifelse(active, "Indirect", "Direct"))
    #
    # Add predicted effect of second intervention
    #
    if (!only.direct) {
      obsdat.i$TE.adj <- ifelse(active, obsdat.i$TE + pred2.i$pred, obsdat.i$TE)
    }
    else {
      obsdat.i$TE.adj <- obsdat.i$TE
      obsdat.i <- obsdat.i[!active, , drop = FALSE]
    }
    #
    # Get the line with predicted values
    #
    seq1 <- seq(0, max(dat.i$.dose1, na.rm = TRUE), length.out = 101)[-1]
    preddat.i <- predict(x, agent1 = i, dose1 = seq1)
    #
    preddat.i <- data.frame(
      agent = i,
      dose1 = preddat.i$dose1,
      pred = preddat.i$pred,
      lower.pred = preddat.i$lower,
      upper.pred = preddat.i$upper
    )
    #
    # Calculation of Benchmark Dose Lower Confidence Limit (BMDL)
    #
    bmdldat.i <- NULL
    #
    if (!is.null(benchmark.threshold)) {
      sel.i <- preddat.i$lower.pred >= benchmark.threshold
      #
      if (any(sel.i))
        bmdldat.i <-
          data.frame(agent = i, bmdl = preddat.i$dose1[which(sel.i)[1]])
    }
    #
    # Calculation of Maximum Effective Dose (MED)
    #
    meddat.i <- NULL
    #
    sel.med.i <- which(abs(diff(preddat.i$pred)) < plateau.threshold)[1]
    if (!is.na(sel.med.i))
      meddat.i <- data.frame(agent = i, med = preddat.i$dose1[sel.med.i + 1])
    #
    dat <- rbind(dat, obsdat.i)
    preddat <- rbind(preddat, preddat.i)
    bmdldat <- rbind(bmdldat, bmdldat.i)
    meddat <- rbind(meddat, meddat.i)
  }
  
  custom_colors <- c("Direct" = col.direct, "Indirect" = col.indirect)
  
  # Create plot
  # 
  p <- ggplot(preddat, aes(x = dose1, y = pred)) +
    geom_ribbon(aes(ymin = lower.pred, ymax = upper.pred),
                fill = "grey40", alpha = 0.4) +
    geom_line(color = col.line) +
    geom_point(data = dat, aes(x = dose1, y = TE.adj, color = type)) +
    scale_color_manual(values = custom_colors) +
    labs(x = "Dose", y = xlab(x$sm, FALSE),
         color = "Observed data comparisons of agent with the reference") +
    theme(plot.title = element_text(size = 10),
          panel.background = element_rect(fill = "white", color = "grey80"),
          panel.grid.major = element_line(color = "grey85", size = 0.3),
          panel.grid.minor = element_line(color = "grey90", size = 0.15),
          plot.background = element_rect(fill = "white", color = NA))
  #
  if (length(coef) > 1)
    p <- p + facet_wrap(~ agent, scales = "free_x")
  #
  if (nrow(bmdldat) > 0)
    p <- p + geom_vline(data = bmdldat, aes(xintercept = bmdl),
                        linetype = "dotted", color = col.bmdl, size = 0.4)
  #
  if (nrow(meddat) > 0)
    p <- p + geom_vline(data = meddat, aes(xintercept = med),
                        linetype = "dotdash", color = col.med, size = 0.4)
  #
  if (legend)
   p <- p + theme(legend.position = "bottom",
                  legend.title = element_text(size = 9),
                  legend.text = element_text(size = 8))
  else
    p <- p + theme(legend.position = "none")
  #
  if (!is.null(ylim))
    p <- p + coord_cartesian(ylim = ylim)
  #
  p
}
