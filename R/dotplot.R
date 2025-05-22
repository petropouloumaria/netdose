#' Dot plot for dose-response data
#'
#' @description
#' Generates a dot plot for dose-response data
#'
#' @param x An object created with \code{\link{netdose}}.
#' @param col The color used for the border of dots.
#' @param fill The color used for the background of dots.
#' @param size A single numeric with the size of the dots.
#' @param drop.reference.group A logical indicating whether to drop the
#'   panel for the reference group.
#' @param ylab A label for the y-axis.
#' @param \dots Additional arguments (ignored).
#'
#' @details
#' The function produces a dot plot of drug doses.
#' 
#' @note The message \emph{'Bin width defaults to 1/30 of the range of the data.
#'   Pick better value with `binwidth`.'} is irrelevant as dot sizes should
#'   be identical. Setting the argument 'binwidth' would result in different
#'   dot sizes for the drugs.
#' 
#' @return No return value.
#' 
#' @keywords hplot
#' 
#' @seealso \code{\link{netdose}}
#'
#' @author Guido Schwarzer <guido.schwarzer@@uniklinik-freiburg.de>
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
#'   append = FALSE)
#' 
#' # DR-NMA with linear dose-response function
#' dr1 <- netdose(TE, seTE, agent1, dose1, agent2,
#'   dose2, studlab,
#'   data = dat)
#'
#' # Dose-response plot
#' dotplot(dr1)
#' 
#' @export dotplot

dotplot <- function(x, col = "black", fill = "orange", size = 2,
                    drop.reference.group = TRUE,
                    ylab = "Dose", ...) {
  # Check class
  chkclass(x, "netdose")
  #
  if (is.null(x$data))
    stop("Data set missing due to using netdose() with ",
         "argument 'keepdata = FALSE'.",
         call. = FALSE)
  #
  chkcolor(col, length = 1)
  chkcolor(fill, length = 1)
  chknumeric(size, min = 0, zero = TRUE, length = 1)
  #
  chklogical(drop.reference.group)
  chkchar(ylab, length = 1)
  
  # Get rid of warnings "no visible binding for global variable"
  #
  studlab <- agent <- .agent <- .agent1 <- .agent2 <-
    dose <- .dose <- .dose1 <- .dose2 <- NULL
  
  # Convert data set to long-arm format
  #
  dat <- x$data %>%
    pivot_longer(cols = c(.agent1, .agent2, .dose1, .dose2),
                 names_to = c(".value", "number"),
                 names_pattern = "(.agent|.dose)([12])") %>%
    select(studlab, .agent, .dose) %>%
    distinct() %>%
    rename(agent = .agent, dose = .dose)
  #
  #
  if (drop.reference.group) {
    dat <- dat %>%
      filter(agent != x$reference.group)
  }
  
  # Create plot
  #
  p <- ggplot(dat, aes(x = 1, y = dose)) +
    geom_dotplot(binaxis = "y",
                 stackdir = "center",
                 method = "histodot",
                 dotsize = size, fill = fill, color = col) +
    facet_wrap(~agent, scales = "free_y") +
    ylab(ylab) +
    theme_minimal() +
    theme(plot.title = element_text(size = 10),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          #
          strip.background = element_rect(fill = "gray80", color = NA),
          #
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  #
  p
}
