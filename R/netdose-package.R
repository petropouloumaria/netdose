#' netdose: Brief overview of network meta-analysis model with dose-response
#' relationships
#'
#' @description
#' R package \bold{netdose} provides methods and graphical tools to conduct
#' the network meta-analysis with dose-response relationships in a
#' frequentist way.
#'
#' @details
#' R package \bold{netdose} is a tool to conduct dose-response network
#' meta-analysis a frequentist way (Petropoulou et al, 2025). The package can
#' implement the dose-response network meta-analysis model (function
#' \code{\link{netdose}}); calculate the predicted values of the dose-response
#' network meta-analysis model (function \code{\link{predict}}); provide
#' dose-response plots (function \code{\link{plot.netdose}})
#' (Petropoulou et al., 2025).
#'
#' Type \code{help(package = "netdose")} for a listing of R functions
#' available in \bold{netdose}.
#'
#' Type \code{citation("netdose")} on how to cite \bold{netdose}
#' in publications.
#'
#' To report problems and bugs, please send an email to Dr. Maria
#' Petropoulou <maria.petropoulou@uniklinik-freiburg.de>.
#'
#' The development version of \bold{netdose} is available on GitHub
#' \url{https://github.com/petropouloumaria/netdose}.
#'
#' @name netdose-package
#'
#' @author Petropoulou Maria <maria.petropoulou@@.uniklinik-freiburg.de>,
#'   Guido Schwarzer <guido.schwarzer@@uniklinik-freiburg.de>
#'
#' @references
#' Petropoulou et al. (2025):
#' Network meta-analysis with dose-response relationships.
#'
#' @keywords package
#'
#' @importFrom meta ci gs
#' @importFrom netmeta netmeta netconnection invmat
#' @importFrom Hmisc rcspline.eval
#' @importFrom MASS ginv
#' @importFrom stats predict quantile median pchisq optimize
#' @importFrom Matrix bdiag
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_ribbon geom_vline
#'   labs theme theme_void element_blank element_text element_rect element_line
#'   annotate scale_color_manual scale_linetype_manual coord_cartesian
#'   ggplotGrob facet_wrap
#' @importFrom grid nullGrob unit unit.c
#' @importFrom gridExtra grid.arrange arrangeGrob
#' @importFrom rlang sym
#' @importFrom dplyr %>% mutate

"_PACKAGE"

NULL
