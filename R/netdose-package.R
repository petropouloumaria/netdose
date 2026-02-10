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
#' meta-analysis a frequentist way (Petropoulou et al, 2026). The package can
#' implement the dose-response network meta-analysis model (function
#' \code{\link{netdose}}); calculate the predicted values of the dose-response
#' network meta-analysis model (function \code{\link{predict}}); provide
#' dose-response plots (function \code{\link{plot.netdose}})
#' (Petropoulou et al., 2026).
#'
#' Type \code{help(package = "netdose")} for a listing of R functions
#' available in \bold{netdose}.
#'
#' Type \code{citation("netdose")} on how to cite \bold{netdose}
#' in publications.
#'
#' To report problems and bugs, please send an email to Dr. Maria
#' Petropoulou <m.petropoulou.a@gmail.com>.
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
#' Petropoulou M, Rücker G, Schwarzer G (2026):
#' Network meta-analysis with dose-response relationships,
#' \emph{BMC Medical Research Methodology},
#' \bold{28}, 17
#'
#' @keywords package
#'
#' @importFrom meta ci gs
#' @importFrom netmeta netmeta netconnection invmat
#' @importFrom Hmisc rcspline.eval
#' @importFrom MASS ginv
#' @importFrom stats predict quantile median pchisq optimize
#' @importFrom Matrix bdiag
#' @importFrom ggplot2 ggplot aes coord_cartesian 
#'   element_blank element_line element_rect element_text
#'   facet_wrap geom_dotplot geom_line geom_point geom_ribbon geom_text
#'   geom_vline guide_axis labs scale_color_manual
#'   scale_x_continuous scale_y_continuous theme theme_minimal
#' @importFrom ggh4x facet_wrap2 facetted_pos_scales
#' @importFrom dplyr %>% arrange count distinct filter group_by if_else
#'   left_join mutate n n_distinct pull rename row_number select
#'   summarise ungroup
#' @importFrom tidyr pivot_longer

"_PACKAGE"

NULL
