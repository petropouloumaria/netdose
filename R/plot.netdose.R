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
#' @param agents Optional character vector specifying which agents to include in the plot. If NULL, all agents will be plotted.
#' @param same.ylim Logical; if TRUE, all plots will have the same y-axis limits. Default is FALSE.
#' @param ylim Optional numeric vector of length 2 specifying the y-axis limits. If NULL, limits are determined automatically (or based on "same.ylim" if TRUE).
#' @param benchmark.threshold Numeric; benchmark response level (e.g., 0.1 for 10 percent). Used to compute Benchmark Dose Lower Confidence Limit (BMDL).
#' @param plateau.threshold Numeric; threshold for identifying the plateau in the dose-response
#'   curve. Defines the minimum absolute change in predicted response between adjacent dose
#'   levels, below which the response is considered stable (i.e., plateau has been reached).
#'   Used to calculate the Maximum Effective Dose (MED). Default: \code{0.0001}.
#' @param ... Additional arguments. Currently ignored, but included for potential future extensions or compatibility with generic plotting functions.
#'
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
#'

plot.netdose <- function(x, pooled = if (x$random) "random" else "common",
                         only.direct = TRUE,
                         col.direct = if (only.direct) "black" else "green",
                         col.indirect = "red",
                         agents = NULL,
                         same.ylim = TRUE,
                         ylim = NULL,
                         benchmark.threshold = 0.01,
                         plateau.threshold = 0.0001,
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
  }
  else {
    coef <- x$TE.drnma.common
  }
  
  # exclude reference group from the plot
  coef <- coef[!names(coef) %in% x$reference.group]
  
  if (!is.null(agents)) {
    coef <- coef[names(coef) %in% agents]
  }
  
  # Create a list to store plots
  plots_list <- list()
  
  all_y_vals <- c()
  
  has_bmdl_or_med <- FALSE
  
  
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
    }
    else {
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
    
    all_y_vals <- c(all_y_vals, linedata$pred, linedata$lower.pred, linedata$upper.pred, dat.obs$TE.adj)
    
    
    benchmark_met <- linedata$pred >= benchmark.threshold
    bmdl <- NA 
    
    if (any(benchmark_met)) {
      bmd_index <- which(benchmark_met)[1]
      bmd_dose <- linedata$dose1[bmd_index]
      
      benchmark_met_l <- linedata$lower.pred >= benchmark.threshold
      
      if (any(benchmark_met_l)) {
        bmd_index_l <- which(benchmark_met_l)[1]
        bmdl <- linedata$dose1[bmd_index_l]
      }
    }
    
    # Calculation of Maximum Effective Dose (MED)
    # 
    delta <- diff(linedata$pred)
    plateau.threshold <- 0.0001
    med_index <- which(abs(delta) < plateau.threshold)[1]
    med <- if (!is.na(med_index)) linedata$dose1[med_index + 1] else NA
    
    
    # Create plot for the current agent
    # 
    plot <- ggplot(linedata, aes(x = dose1, y = pred)) +
      geom_ribbon(aes(ymin = lower.pred, ymax = upper.pred),
                  fill = "grey40", alpha = 0.4) +
      geom_line(color = "blue") +
      labs(y =  NULL, x = NULL, title = i,
           color = "Observed data comparisons of agent with the reference") +
      theme(plot.title = element_text(size = 10),
            panel.background = element_rect(fill = "white", color = "grey80"),
            panel.grid.major = element_line(color = "grey85", size = 0.3),
            panel.grid.minor = element_line(color = "grey90", size = 0.15),
            plot.background = element_rect(fill = "white", color = NA),
            legend.position = "bottom",
            legend.title = element_text(size = 9),
            legend.text = element_text(size = 8))
    
    # Add BMDL to the plot    
    if (any(benchmark_met == TRUE)  && !is.na(bmdl))  {
      has_bmdl_or_med <- TRUE
      plot <- plot +
        geom_vline(xintercept = bmdl, linetype = "dotted",
                   color = "purple", size = 0.4) +
        annotate("text", x = bmdl, y = linedata$lower.pred[bmd_index] - 0.001,
                 label = paste0(round(bmdl, 3)), angle = 0, size = 2.5,
                 color = "purple")
    }
    
    # Add MED to the plot
    if (!is.na(med) && !is.na(bmdl) && med >= bmdl) {
      has_bmdl_or_med <- TRUE
      plot <- plot +
        geom_vline(xintercept = med, linetype = "dotdash", color = "gray40", size = 0.4) +
        annotate("text",
                 x = med,
                 y = linedata$upper.pred[med_index] + 0.001,
                 label = paste0(round(med, 3)),
                 angle = 0, size = 2.5, color = "gray40")
    }
    
    if (nrow(dat.obs) > 0) {
      if (only.direct) {
        dat.obs$Comparison <- "Direct"
        plot <- plot +
          geom_point(data = dat.obs,
                     aes(x = dose1, y = TE.adj,
                         color = !!sym("Comparison")), size = 1) +
          scale_color_manual(values = c("Direct" = col.direct))
      }
      else {
        dat.obs$Comparison <- ifelse(active, "Indirect", "Direct")
        plot <- plot +
          geom_point(data = dat.obs,
                     aes(x = dose1, y = TE.adj,
                         color = !!sym("Comparison")), size = 1) +
          scale_color_manual(values = c("Direct" = col.direct,
                                        "Indirect" = col.indirect))
      }
    }
    plots_list[[i]] <- plot
  }
  
  
  label = c("BMDL", "MED")
  # Create dummy data for legend (BMDL and MED)
  if (has_bmdl_or_med) {
    vline_legend_df <- data.frame(
      x = c(1, 2),
      label = c("BMDL", "MED")
    )
    
    # Dummy plot to extract legend
    legend_bmd_med <- ggplot() +
      geom_vline(data = vline_legend_df,
                 aes(xintercept = x, color = label, linetype = label),
                 size = 0.6) +
      scale_color_manual(values = c("BMDL" = "purple", "MED" = "gray40")) +
      scale_linetype_manual(values = c("BMDL" = "dotted", "MED" = "dotdash")) +
      theme_void() +
      theme(legend.position = "bottom",
            legend.title = element_blank(),
            legend.text = element_text(size = 9))
    
    # Extract only the legend grob
    g_legend <- ggplotGrob(legend_bmd_med)
    legend_index <-
      which(sapply(g_legend$grobs, function(x) x$name) == "guide-box")
    
    if (length(legend_index) > 0) {
      shared_bmd_med_legend <- g_legend$grobs[[legend_index]]
    }
    else {
      shared_bmd_med_legend <- nullGrob()
    }
  }
  else {
    shared_bmd_med_legend <- nullGrob()
  }
  
  
  if (!is.null(ylim)) {
    plots_list <-
      lapply(plots_list, function(p) p + coord_cartesian(ylim = ylim))
  }
  else if (same.ylim) {
    ylim_global <- range(all_y_vals, na.rm = TRUE)
    plots_list <-
      lapply(plots_list, function(p) p + coord_cartesian(ylim = ylim_global))
  }
  
  if (x$method == "fp1") {
    if (is.null(x$param)) {
      x$param <- -0.5
    }
    method <- paste0("FP1 (p=", x$param, ")")
  }
  else if (x$method == "fp2") {
    if (is.null(x$param)) {
      x$param <- c(-0.5, 1)
    }
    method <- paste0("FP2 (p1=", x$param[1], ", p2=", x$param[2], ")")
  }
  else if (x$method == "rcs") {
    method <- paste0("RCS (p=", paste(x$param, collapse = ", "), ")")
  }
  else {
    method <- paste0(toupper(substring(x$method, 1, 1)), substring(x$method, 2))
  }
  
  
  legend_plot <- plots_list[[1]] + theme(legend.position = "bottom")
  
  g <- ggplotGrob(legend_plot)
  legend_index <- which(sapply(g$grobs, function(x) x$name) == "guide-box")
  
  if (length(legend_index) > 0) {
    shared_legend <- g$grobs[[legend_index]]
  }
  else {
    shared_legend <- nullGrob() # Empty placeholder grob
  }
  
  plots_nolegend <-
    lapply(plots_list, function(p) p + theme(legend.position = "none"))
  
  plots_grid <- arrangeGrob(
    grobs = plots_nolegend,
    ncol = min(3, length(plots_nolegend)),
    top = paste(method, "dose-response plot"),
    left = xlab(x$sm, FALSE),
    bottom = "Dose"
  )
  
  grid.arrange(
    plots_grid,
    shared_legend,
    shared_bmd_med_legend,
    ncol = 1,
    heights = unit.c(
      unit(1, "npc") - unit(3, "lines"),
      unit(1.5, "lines"),
      unit(1.5, "lines")
    )
  )
}
