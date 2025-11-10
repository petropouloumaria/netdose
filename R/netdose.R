#' Network meta-analysis with dose-response relationships
#'
#' @description
#' The `netdose` function performs a dose-response network meta-analysis in a
#' frequentist way. It accepts a dataset with study-level data, constructs a
#' design matrix for the dose-response model, and computes treatment effects
#' under the common and random effects models. The function supports multiple
#' dose-response relationship modelling approaches, including linear,
#' exponential, quadratic, restricted cubic splines (rcs), and fractional
#' polynomials (fp1, fp2).
#'
#' @param TE Estimate of treatment effect, i.e. difference between
#'   first and second treatment (e.g. log odds ratio, mean difference,
#'   or log hazard ratio). Or an R object created with
#'   \code{\link[meta]{pairwise}}.
#' @param seTE Standard error of treatment estimate.
#' @param agent1 Agents corresponding to the first treatment in each comparison.
#' @param dose1 Doses for the first treatment in each comparison.
#' @param agent2 Agents corresponding to the second treatment in each
#'   comparison.
#' @param dose2 Doses for the second treatment in each comparison.
#' @param studlab An optional - but important! - vector with study labels.
#' @param data An optional data frame containing the study information.
#' @param subset An optional vector specifying a subset of studies to
#'   be used. The default is `NULL`.
#' @param n1 Numeric. Optional. Sample sizes for the first treatment in
#'   each comparison.
#' @param n2 Numeric. Optional. Sample sizes for the second treatment in
#'   each comparison.
#' @param event1 Numeric. Optional. Number of events for the first treatment
#'   in each comparison.
#' @param event2 Numeric. Optional. Number of events for the second treatment
#'   in each comparison.
#' @param sm  A character string indicating underlying summary measure,
#'   e.g., \code{"RD"}, \code{"RR"}, \code{"OR"}, \code{"ASD"},
#'   \code{"HR"}, \code{"MD"}, \code{"SMD"}, or \code{"ROM"}.
#' @param common A logical indicating whether a common effects dose-response
#'   network meta-analysis should be conducted. The default is \code{TRUE}.
#' @param random A logical indicating whether a random effects dose-response
#'   network meta-analysis should be conducted. The default is \code{TRUE}.
#' @param tau An optional value for the square-root of the
#'   between-study variance \eqn{\tau^2}.
#' @param method An optional character string specifying the method to
#'  be used for the dose-response relationship. Either, "linear", "exponential",
#'  "quadratic", "rcs", "fp1", or "fp2", can be abbreviated (see Details).
#' @param param A numeric vector specifying the parameters for some
#'   dose-response functions (see Details).
#' @param reference.group Reference agent (first agent with dose 0 is used
#'   if argument is missing).
#' @param common.dose A named vector with the common dose for each agent
#'   in the network (see Examples). The median dose is used for each agent if
#'   this argument is not provided.
#' @param level The level used to calculate confidence intervals for
#' individual comparisons.
#' @param backtransf A logical indicating whether results should be
#' back transformed in printouts and forest plots. If
#' \code{backtransf = TRUE}, results for \code{sm = "OR"} are
#' presented as odds ratios rather than log odds ratios, for
#' example.
#' @param tol.multiarm A numeric for the tolerance for consistency of
#'   treatment estimates in multi-arm studies which are consistent by
#'   design (only considered for standard network meta-analysis model).
#' @param tol.multiarm.se A numeric for the tolerance for consistency
#'   of standard errors in multi-arm studies which are consistent by
#'   design (only considered for standard network meta-analysis model).
#' @param details.chkmultiarm A logical indicating whether treatment
#'   estimates and / or variances of multi-arm studies with
#'   inconsistent results or negative multi-arm variances should be
#'   printed (only considered for standard network meta-analysis model).
#' @param keepdata A logical indicating whether original data(set)
#'   should be kept in netdose object.
#' @param warn A logical indicating whether warnings should be printed
#'   (e.g., if studies are excluded from network meta-analysis due to zero
#'   standard errors).
#' @param func.inverse R function used to calculate the pseudoinverse
#'   of the Laplacian matrix L.
#'
#' @details
#' The dose-response network meta-analysis (DR-NMA) has been implemented
#' by modelling different dose-response functions, as described by
#' Mandema et al. 2005 and Mawdsley et al. 2016 and by using restricted cubic
#' splines (Hamza et al. 2024) in a Bayesian setting.
#'
#' The function \code{netdose} conducts a dose-response network meta-analysis
#' with a variety of dose-response functions (such as the linear, exponential,
#' fractional polynomials and restricted cubic splines) in a frequentist way
#' as described in Petropoulou et al. (2025).
#'
#' The following dose-response functions are available:
#' \itemize{
#' \item Linear dose-response relationship (\code{method = "linear"})
#' \item Exponential dose-response relationship (\code{method = "exponential"})
#' \item Quadratic polynomial dose-response relationship
#'   (\code{method = "quadratic"})
#' \item Restricted cubic splines (\code{method = "rcs"})
#' \item Fractional polynomial (order 1) (\code{method = "fp1"})
#' \item Fractional polynomial (order 2) (\code{method = "fp2"})
#' }
#' By default, a linear dose-response relationship is assumed.
#'
#' The parameters for the selected dose-response function can be specified
#' using argument \code{param}: a numeric vector specifying the percentiles to
#' set the knots for the restricted cubic splines (default: knots at the 10th,
#' 50th, and 90th percentile), a single numeric specifying the power of the
#' fractional polynomial with order 1 (default: -0.5), or a numeric vector of
#' length 2 specifying the first and second power of a fractional polynomial
#' with order 2 (default: -0.5 and -0.5). The input for argument \code{param}
#' is ignored for a linear, exponential or quadratic polynomial dose-response
#' relationship.
#'
#' @return
#' An object of class \code{netdose}; a list containing the
#' following components:
#' \item{studlab}{Study labels.}
#' \item{agent1}{Label/Agents corresponding to the first treatment in each
#'   comparison.}
#' \item{agent2}{Label/Agents corresponding to the second treatment in each
#'   comparison.}
#' \item{dose1}{Doses for the first treatment in each comparison.}
#' \item{dose2}{Doses for the second treatment in each comparison.}
#' \item{treat1}{Label/First treatment in each comparison.}
#' \item{treat2}{Label/Second treatment in each comparison.}
#'
#' \item{TE}{Estimate of treatment effect, i.e. difference between
#'   first and second treatment.}
#' \item{seTE}{Standard error of treatment estimate.}
#' \item{seTE.adj.common, seTE.adj.random}{Standard error of treatment
#'   estimate, adjusted for multi-arm studies.}
#'
#' \item{k}{Total number of studies.}
#' \item{m}{Total number of pairwise comparisons.}
#' \item{a}{Total number of agents.}
#' \item{n}{Total number of treatments.}
#' \item{trts}{Treatments included in the dataset in alphabetic order.}
#' \item{agents}{Agents included in dose-response network meta-analysis in
#' alphabetic order.}
#' \item{inactive}{Identifier for the reference group or inactive treatment.}
#'
#' \item{common.dose}{Common dose value used in the analysis, if specified.}
#'
#' **Common/Random effects model results:**
#'  \item{TE.common, TE.random}{Matrix with overall
#'   treatment effects estimated by the dose-response (common and random
#'   effects) model.}
#' \item{seTE.common, seTE.random}{Matrix with
#'   standard errors estimated by the dose-response (common and random
#'   effects) model.}
#' \item{lower.common, upper.common, lower.random,
#'   upper.random}{Matrices with lower and upper
#'   confidence interval limits estimated by the dose-response (common and
#'   random effects) model.}
#' \item{statistic.common, pval.common, statistic.random,
#'   pval.random}{Matrices with z-values and
#'   p-values for test of overall effect estimated by the dose-response
#'   (common and random effects) model.}
#' \item{TE.drnma.common}{A vector of dose-response effects (common
#'   and random effects model).}
#' \item{seTE.drnma.common, seTE.drnma.random}{A vector with corresponding
#'   standard errors (common and random effects model).}
#' \item{lower.drnma.common, lower.drnma.random}{A vector with lower
#'   confidence limits for dose-response treatment estimates (common and
#'   random effects model).}
#' \item{upper.drnma.common, upper.drnma.random}{A vector with upper
#'  confidence limits for dose-response treatment estimates (common and
#'  random effects model).}
#' \item{statistic.drnma.common, statistic.drnma.random}{A vector with
#'   z-values for the overall dose-response effects (common and random
#'   effects model).}
#' \item{pval.drnma.common, pval.drnma.random}{A vector with p-values for
#'   the overall dose-response effects (common and random effects
#'   model).}
#'
#' **Heterogeneity and goodness-of-fit statistics:**
#' \item{Q}{Overall heterogeneity / inconsistency statistic for
#'   dose-response network meta-analysis.}
#' \item{df.Q}{Degrees of freedom for test of heterogeneity /
#'   inconsistency for dose-response network meta-analysis.}
#' \item{pval.Q}{P-value for test of heterogeneity /
#'   inconsistency for dose-response network meta-analysis.}
#' \item{tau}{Square-root of between-study variance with DerSimonian
#'   and Laird method for dose-response network meta-analysis.}
#' \item{tauml}{Square-root of between-study variance with Maximum
#'   likelihood method for dose-response network meta-analysis.}
#' \item{I2, lower.I2, upper.I2}{I-squared, lower and upper confidence
#'   limits.}
#' \item{Q.lump}{Overall heterogeneity / inconsistency statistic
#'   (NMA, lumping approach).}
#' \item{df.Q.lump}{Degrees of freedom for test of heterogeneity /
#'   inconsistency (NMA, lumping approach).}
#' \item{pval.Q.lump}{P-value for test of heterogeneity /
#'   inconsistency (NMA, lumping approach).}
#' \item{Q.split}{Overall heterogeneity / inconsistency statistic
#'   (NMA, splitting approach).}
#' \item{df.Q.split}{Degrees of freedom for test of heterogeneity /
#'   inconsistency (NMA, splitting approach).}
#' \item{pval.Q.split}{P-value for test of heterogeneity /
#'   inconsistency (NMA, splitting approach).}
#'
#' \item{B.matrix}{Edge-vertex incidence matrix.}
#' \item{D_obs.matrix}{Matrix with observed doses.}
#' \item{D.matrix}{Matrix with transformed doses.}
#' \item{X.matrix}{Design matrix for dose-response network meta-analysis.}
#'
#' \item{sm}{Summary measure used in the analysis.}
#' \item{level}{Level used to calculate confidence intervals for
#'   individual comparisons.}
#' \item{common}{A logical indicating whether a common effects dose-response
#'   network meta-analysis should be conducted.}
#' \item{random}{A logical indicating whether a random effects dose-response
#'   network meta-analysis should be conducted.}
#' \item{method}{Method used for the dose-response relationship.}
#'
#' \item{reference.group}{Reference agent.}
#' \item{Q.to.df.ratio}{Q to df ratio, i.e, Q/df.Q.}
#' \item{func.inverse}{Function used to calculate the pseudoinverse of
#'   the Laplacian matrix L.}
#' \item{backtransf}{A logical indicating whether results should be
#'   back transformed in printouts and forest plots.}
#'
#' \item{data}{Data frame containing the study information.}
#'
#' @author Maria Petropoulou <m.petropoulou.a@gmail.com>,
#'   Guido Schwarzer <guido.schwarzer@@uniklinik-freiburg.de>
#'
#' @references
#' Mandema JW, Cox EJ (2005):
#' Therapeutic benefit of eletriptan compared to sumatriptan for the acute
#' relief of migraine pain--results of a model-based meta-analysis that accounts
#' for encapsulation.
#' \emph{Cephalalgia},
#' \bold{25}, 715--25
#'
#' Mawdsley D, Bennetts M, Dias S, Boucher M, Welton N (2016):
#' Model-Based Network Meta-Analysis: A Framework for Evidence Synthesis of
#' Clinical Trial Data.
#' \emph{PT Pharmacometrics & Systems Pharmacology},
#' \bold{5}, 393--401
#'
#' Hamza T, Furukawa TA, Orsin N, Cipriani A, Iglesias CP, Salanti G (2024):
#' A dose-effect network meta-analysis model with application in antidepressants
#' using restricted cubic splines.
#' \emph{Statistical Methods in Medical Research},
#' \bold{33}, 1461--72
#'
#' Petropoulou et al. (2025):
#' Network meta-analysis with dose-response relationships.
#'
#' @examples
#' # Use a subset of 3 studies from anesthesia data
#' anesthesia_subset <-
#'   subset(anesthesia, study %in% unique(anesthesia$study)[1:3])
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
#' # Perform DR-NMA with a linear dose-response function
#' dr1 <- netdose(
#'   TE, seTE, agent1, dose1, agent2,
#'   dose2, studlab,
#'   data = dat
#' )
#' 
#' \donttest{
#' # DR-NMA with FP1 dose-response function with p = -0.5
#' dr_fp1 <- netdose(TE, seTE, agent1, dose1, agent2,
#'   dose2, studlab,
#'   data = dat,
#'   method = "fp1"
#' )
#'
#' # DR-NMA with FP1 dose-response function with p = 0.5
#' dr_fp1_p0.5 <- netdose(TE, seTE, agent1, dose1, agent2,
#'   dose2, studlab,
#'   data = dat,
#'   method = "fp1", param = 0.5
#' )
#'
#' # DR-NMA with RCS dose-response function with knots at 10th, 50th and 90th percentiles
#' dr_rcs <- netdose(TE, seTE, agent1, dose1, agent2,
#'   dose2, studlab,
#'   data = dat,
#'   method = "rcs"
#' )
#'
#' # DR-NMA with RCS dose-response function with knots at 25th, 50th and 100th percentiles
#' dr_rcs2 <- netdose(TE, seTE, agent1, dose1, agent2,
#'   dose2, studlab,
#'   data = dat,
#'   method = "rcs", p = c(0.25, 0.50, 1),
#' )
#' }
#'
#' @export netdose

netdose <- function(TE, seTE, agent1, dose1, agent2, dose2, studlab,
                    data = NULL, subset = NULL,
                    #
                    n1 = NULL, n2 = NULL,
                    event1 = NULL, event2 = NULL,
                    #
                    sm,
                    common = gs("common"),
                    random = gs("random") | !is.null(tau),
                    tau = NULL,
                    #
                    method = "linear",
                    param = NULL,
                    #
                    reference.group,
                    #
                    common.dose = NULL,
                    #
                    level = gs("level.comb"),
                    backtransf = gs("backtransf"),
                    #
                    tol.multiarm = 0.001,
                    tol.multiarm.se = NULL,
                    details.chkmultiarm = FALSE,
                    #
                    func.inverse = invmat,
                    #
                    keepdata = gs("keepdata"),
                    #
                    warn = TRUE) {
  #
  #
  # (1) Check arguments
  #
  #
  chklogical(common)
  chklogical(random)
  #
  if (!is.null(tau)) {
    chknumeric(tau, min = 0, length = 1)
  }
  #
  method <-
    setchar(method, c(
      "linear", "exponential", "quadratic", "rcs",
      "fp1", "fp2"
    ))
  #
  if (method == "rcs")
    chknumeric(param)
  else if (method == "fp2")
    chknumeric(param, length = 2)
  else
    chknumeric(param, length = 1)
  #
  chklevel(level)
  chklogical(backtransf)
  #
  chknumeric(tol.multiarm, min = 0, length = 1)
  if (!is.null(tol.multiarm.se)) {
    chknumeric(tol.multiarm.se, min = 0, length = 1)
  }
  chklogical(details.chkmultiarm)
  #
  chklogical(keepdata)
  chklogical(warn)
  
  
  #
  #
  # (2) Read data
  #
  #
  
  nulldata <- is.null(data)
  sfsp <- sys.frame(sys.parent())
  mc <- match.call()
  #
  TE <- catch("TE", mc, data, sfsp)
  #
  missing.reference.group <- missing(reference.group)
  missing.reference.group.pairwise <- FALSE
  #
  if (is.data.frame(TE) & !is.null(attr(TE, "pairwise"))) {
    is.pairwise <- TRUE
    #
    sm <- attr(TE, "sm")
    if (missing.reference.group) {
      missing.reference.group.pairwise <- TRUE
      reference.group <- attr(TE, "reference.group")
      #
      if (is.null(reference.group))
        reference.group <- ""
      else
        missing.reference.group <- FALSE
    }
    #
    keep.all.comparisons <- attr(TE, "keep.all.comparisons")
    if (!is.null(keep.all.comparisons) && !keep.all.comparisons) {
      stop("First argument is a pairwise object created with ",
           "'keep.all.comparisons = FALSE'.",
           call. = TRUE
      )
    }
    #
    seTE <- TE$seTE
    studlab <- TE$studlab
    #
    if (!is.null(TE$n1)) {
      n1 <- TE$n1
    }
    if (!is.null(TE$n2)) {
      n2 <- TE$n2
    }
    if (!is.null(TE$event1)) {
      event1 <- TE$event1
    }
    if (!is.null(TE$event2)) {
      event2 <- TE$event2
    }
    #
    pairdata <- TE
    data <- TE
    #
    agent1 <- catch("agent1", mc, data, sfsp)
    dose1 <- catch("dose1", mc, data, sfsp)
    #
    agent2 <- catch("agent2", mc, data, sfsp)
    dose2 <- catch("dose2", mc, data, sfsp)
    #
    TE <- TE$TE
  }
  else {
    is.pairwise <- FALSE
    if (missing(sm)) {
      if (!is.null(data) && !is.null(attr(data, "sm")))
        sm <- attr(data, "sm")
      else
        sm <- ""
    }
    #
    seTE <- catch("seTE", mc, data, sfsp)
    #
    agent1 <- catch("agent1", mc, data, sfsp)
    dose1 <- catch("dose1", mc, data, sfsp)
    #
    agent2 <- catch("agent2", mc, data, sfsp)
    dose2 <- catch("dose2", mc, data, sfsp)
    #
    studlab <- catch("studlab", mc, data, sfsp)
    #
    n1 <- catch("n1", mc, data, sfsp)
    n2 <- catch("n2", mc, data, sfsp)
    #
    event1 <- catch("event1", mc, data, sfsp)
    event2 <- catch("event2", mc, data, sfsp)
  }
  #
  chknumeric(TE)
  chknumeric(seTE)
  #
  if (!any(!is.na(TE) & !is.na(seTE))) {
    stop("Missing data for estimates (argument 'TE') and ",
         "standard errors (argument 'seTE') in all studies.\n  ",
         "No network meta-analysis possible.",
         call. = FALSE
    )
  }
  #
  k.Comp <- length(TE)
  #
  chklength(seTE, k.Comp, "TE")
  chklength(agent1, k.Comp, "TE")
  chklength(dose1, k.Comp, "TE")
  chklength(agent2, k.Comp, "TE")
  chklength(dose2, k.Comp, "TE")
  #
  if (is.null(studlab))
    studlab <- seq_along(TE)
  else
    chklength(studlab, k.Comp, "TE")
  #
  if (is.factor(agent1))
    agent1 <- as.character(agent1)
  #
  if (is.factor(agent2))
    agent2 <- as.character(agent2)
  #
  # Remove leading and trailing whitespace
  #
  agent1 <- rmSpace(rmSpace(agent1, end = TRUE))
  agent2 <- rmSpace(rmSpace(agent2, end = TRUE))
  #
  if (length(studlab) == 0) {
    if (warn) {
      warning("No information given for argument 'studlab'. ",
              "Assuming that comparisons are from independent studies.",
              call. = FALSE
      )
    }
    studlab <- seq(along = TE)
  }
  studlab <- as.character(studlab)
  #
  subset <- catch("subset", mc, data, sfsp)
  missing.subset <- is.null(subset)
  #
  if (!is.null(event1) & !is.null(event2))
    available.events <- TRUE
  else
    available.events <- FALSE
  #
  if (!is.null(n1) & !is.null(n2))
    available.n <- TRUE
  else
    available.n <- FALSE
  #
  treat1 <- paste(agent1, dose1)
  treat2 <- paste(agent2, dose2)
  
  
  #
  #
  # (2b) Store complete dataset in list object data
  #      (if argument keepdata is TRUE)
  #
  #
  if (keepdata) {
    if (nulldata & !is.pairwise) {
      data <- data.frame(.studlab = studlab, stringsAsFactors = FALSE)
    }
    else if (nulldata & is.pairwise) {
      data <- pairdata
      data$.studlab <- studlab
    }
    else
      data$.studlab <- studlab
    #
    data$.order <- seq_along(studlab)
    #
    data$.agent1 <- agent1
    data$.dose1 <- dose1
    data$.treat1 <- treat1
    #
    data$.agent2 <- agent2
    data$.dose2 <- dose2
    data$.treat2 <- treat2
    #
    data$.TE <- TE
    data$.seTE <- seTE
    #
    data$.event1 <- event1
    data$.n1 <- n1
    data$.event2 <- event2
    data$.n2 <- n2
    #
    # Check for correct agent order within comparison
    #
    wo <- data$.agent1 > data$.agent2
    #
    if (any(wo)) {
      data$.TE[wo] <- -data$.TE[wo]
      #
      tagent1 <- data$.agent1
      data$.agent1[wo] <- data$.agent2[wo]
      data$.agent2[wo] <- tagent1[wo]
      #
      tdose1 <- data$.dose1
      data$.dose1[wo] <- data$.dose2[wo]
      data$.dose2[wo] <- tdose1[wo]
      #
      ttreat1 <- data$.treat1
      data$.treat1[wo] <- data$.treat2[wo]
      data$.treat2[wo] <- ttreat1[wo]
      #
      if (isCol(data, ".n1") & isCol(data, ".n2")) {
        tn1 <- data$.n1
        data$.n1[wo] <- data$.n2[wo]
        data$.n2[wo] <- tn1[wo]
      }
      #
      if (isCol(data, ".event1") & isCol(data, ".event2")) {
        tevent1 <- data$.event1
        data$.event1[wo] <- data$.event2[wo]
        data$.event2[wo] <- tevent1[wo]
      }
    }
    #
    if (!missing.subset) {
      if (length(subset) == dim(data)[1]) {
        data$.subset <- subset
      }
      else {
        data$.subset <- FALSE
        data$.subset[subset] <- TRUE
      }
    }
  }
  
  
  #
  #
  # (3) Use subset for analysis
  #
  #
  if (!missing.subset) {
    if ((is.logical(subset) & (sum(subset) > k.Comp)) ||
        (length(subset) > k.Comp)) {
      stop("Length of subset is larger than number of studies.",
           call. = FALSE
      )
    }
    #
    TE <- TE[subset]
    seTE <- seTE[subset]
    #
    agent1 <- agent1[subset]
    agent2 <- agent2[subset]
    #
    dose1 <- dose1[subset]
    dose2 <- dose2[subset]
    #
    treat1 <- treat1[subset]
    treat2 <- treat2[subset]
    #
    studlab <- studlab[subset]
    #
    if (!is.null(n1)) {
      n1 <- n1[subset]
    }
    if (!is.null(n2)) {
      n2 <- n2[subset]
    }
    if (!is.null(event1)) {
      event1 <- event1[subset]
    }
    if (!is.null(event2)) {
      event2 <- event2[subset]
    }
  }
  #
  trts <- sort(unique(c(treat1, treat2)))
  n <- length(trts)
  #
  agents <- sort(unique(c(agent1, agent2)))
  #
  if (!is.null(common.dose)) {
    chklength(common.dose, length(agents),
              text =
                paste(
                  "Argument 'common.dose' must be a named vector of",
                  "length", length(agents),
                  "(number of agents in network)."
                )
    )
    #
    names(common.dose) <- setchar(names(common.dose), agents)
  }
  else {
    by_dose <- by(c(dose1, dose2), c(agent1, agent2), median, na.rm = TRUE)
    #
    common.dose <- as.vector(by_dose)
    names(common.dose) <- names(by_dose)
  }
  
  
  #
  #
  # (4) Additional checks
  #
  #
  
  if (any(agent1 == agent2 & dose1 == dose2)) {
    stop("Same combination of agent and dose in at least one comparison.",
         call. = FALSE
    )
  }
  #
  # Check NAs in estimates or standard errors and zero standard errors
  #
  excl <- is.na(TE) | is.na(seTE) | seTE <= 0
  #
  if (any(excl)) {
    if (keepdata) {
      data$.excl <- excl
    }
    #
    dat.NAs <- data.frame(
      studlab = studlab[excl],
      agent1 = agent1[excl],
      dose1 = dose1[excl],
      agent2 = agent2[excl],
      dose2 = dose2[excl],
      TE = format(round(TE[excl], 4)),
      seTE = format(round(seTE[excl], 4)),
      stringsAsFactors = FALSE
    )
    if (warn) {
      warning("Comparison",
              if (sum(excl) > 1) "s",
              " with missing TE / seTE or zero seTE not considered ",
              "in network meta-analysis.",
              call. = FALSE
      )
    }
    if (warn) {
      cat("Comparison",
          if (sum(excl) > 1) "s",
          " not considered in network meta-analysis:\n",
          sep = ""
      )
      prmatrix(dat.NAs,
               quote = FALSE, right = TRUE,
               rowlab = rep("", sum(excl))
      )
      cat("\n")
    }
    #
    studlab <- studlab[!excl]
    agent1 <- agent1[!excl]
    agent2 <- agent2[!excl]
    dose1 <- dose1[!excl]
    dose2 <- dose2[!excl]
    treat1 <- treat1[!excl]
    treat2 <- treat2[!excl]
    TE <- TE[!excl]
    seTE <- seTE[!excl]
    #
    if (!is.null(n1)) {
      n1 <- n1[!excl]
    }
    if (!is.null(n2)) {
      n2 <- n2[!excl]
    }
    if (!is.null(event1)) {
      event1 <- event1[!excl]
    }
    if (!is.null(event2)) {
      event2 <- event2[!excl]
    }
  }
  #
  # Check NAs in dose1 or dose2
  #
  excl <- is.na(dose1) | is.na(dose2)
  #
  if (any(excl)) {
    if (keepdata) {
      data$.excl <- excl
    }
    #
    dat.NAs <- data.frame(
      studlab = studlab[excl],
      agent1 = agent1[excl],
      dose1 = dose1[excl],
      agent2 = agent2[excl],
      dose2 = dose2[excl],
      TE = format(round(TE[excl], 4)),
      seTE = format(round(seTE[excl], 4)),
      stringsAsFactors = FALSE
    )
    if (warn) {
      warning("Comparison",
              if (sum(excl) > 1) "s",
              " with missing dose information not considered ",
              "in network meta-analysis.",
              call. = FALSE
      )
    }
    if (warn) {
      cat("Comparison",
          if (sum(excl) > 1) "s",
          " not considered in network meta-analysis:\n",
          sep = ""
      )
      prmatrix(dat.NAs,
               quote = FALSE, right = TRUE,
               rowlab = rep("", sum(excl))
      )
      cat("\n")
    }
    #
    studlab <- studlab[!excl]
    agent1 <- agent1[!excl]
    agent2 <- agent2[!excl]
    dose1 <- dose1[!excl]
    dose2 <- dose2[!excl]
    treat1 <- treat1[!excl]
    treat2 <- treat2[!excl]
    TE <- TE[!excl]
    seTE <- seTE[!excl]
    #
    if (!is.null(n1)) {
      n1 <- n1[!excl]
    }
    if (!is.null(n2)) {
      n2 <- n2[!excl]
    }
    if (!is.null(event1)) {
      event1 <- event1[!excl]
    }
    if (!is.null(event2)) {
      event2 <- event2[!excl]
    }
  }
  #
  # Check for correct agent order within comparison
  #
  wo <- agent1 > agent2
  #
  if (any(wo)) {
    TE[wo] <- -TE[wo]
    #
    tagent1 <- agent1
    agent1[wo] <- agent2[wo]
    agent2[wo] <- tagent1[wo]
    #
    tdose1 <- dose1
    dose1[wo] <- dose2[wo]
    dose2[wo] <- tdose1[wo]
    #
    ttreat1 <- treat1
    treat1[wo] <- treat2[wo]
    treat2[wo] <- ttreat1[wo]
    #
    if (available.n) {
      tn1 <- n1
      n1[wo] <- n2[wo]
      n2[wo] <- tn1[wo]
    }
    #
    if (available.events) {
      tevent1 <- event1
      event1[wo] <- event2[wo]
      event2[wo] <- tevent1[wo]
    }
  }
  #
  # Inactive agents, i.e., dose = 0
  #
  by_zero <- by(c(dose1, dose2) == 0, c(agent1, agent2), all)
  inactive <- names(by_zero)[as.vector(by_zero)]
  #
  # Check for agents with dose = 0 and dose > 0
  #
  by_withzero <- by(c(dose1, dose2) == 0, c(agent1, agent2), any)
  withzero <- names(by_withzero)[as.vector(by_withzero)]
  #
  active_and_inactive <- withzero[!(withzero %in% inactive)]
  n.ai <- length(active_and_inactive)
  #
  if (n.ai > 0) {
    stop("The following agent",
         if (n.ai > 1) "s have " else " has ",
         "doses equal to zero and larger than zero: ",
         paste0("'", active_and_inactive, "'", collapse = ", "),
         call. = FALSE
    )
  }
  #
  # Check value for reference group
  #
  labels <- sort(unique(agents))
  #
  if (missing.reference.group | missing.reference.group.pairwise) {
    if (length(inactive) > 0) {
      reference.group <- inactive[1]
    }
    else {
      go.on <- TRUE
      i <- 0
      while (go.on) {
        i <- i + 1
        sel.i <-
          !is.na(TE) & !is.na(seTE) &
          (agent1 == labels[i] | agent2 == labels[i])
        if (sum(sel.i) > 0) {
          go.on <- FALSE
          reference.group <- labels[i]
        }
        else if (i == length(labels)) {
          go.on <- FALSE
          reference.group <- ""
        }
      }
    }
  }
  #
  if (reference.group != "") {
    reference.group <- setref(reference.group, labels)
  }
  
  
  #
  #
  # (5) Create design matrix X
  #
  #
  
  if (method == "linear") {
    param <- NULL
    #
    M <- createXd1(agent1, dose1, agent2, dose2, studlab, g = dose2dose)
  }
  else if (method == "exponential") {
    if (is.null(param)) {
      param <- 1
    }
    #
    M <- createXd1(agent1, dose1, agent2, dose2, studlab,
                   g = dose2exp, param = param)
  }
  else if (method == "fp1") { # fractional polynomial (order 1)
    if (is.null(param)) {
      param <- -0.5
    }
    #
    M <- createXd1(agent1, dose1, agent2, dose2, studlab,
                   g = dose2fp, param = param
    )
  }
  else if (method == "quadratic") { # quadratic polynomial
    if (is.null(param)) {
       param <- c(1, 2)
    }
    #
    M <- createXd2(agent1, dose1, agent2, dose2, studlab,
                   g1 = dose2dose, g2 = dose2poly, param = param
    )
  }
  else if (method == "fp2") { # fractional polynomial (order 2)
    if (is.null(param)) {
      param <- c(-0.5, 1)
    }
    #
    M <- createXd2(agent1, dose1, agent2, dose2, studlab,
                   g1 = dose2fp, g2 = dose2fp, param = param
    )
  }
  else if (method == "rcs") { # restricted cubic splines
    if (is.null(param)) { # Harrell's suggestion (fixed sample quantiles)
      param <- c(0.10, 0.50, 0.90)
    }
    # Alternative choice: c(0.25, 0.50, 1)
    #
    M <- createXd_rcs(agent1, dose1, agent2, dose2, studlab, param = param)
  }
  #
  B.matrix <- M$B.matrix
  D_obs.matrix <- M$D_obs.matrix
  D.matrix <- M$D.matrix
  X.matrix <- M$X.matrix
  
  
  #
  #
  # (6) Conduct DR-NMA
  #
  #
  
  p0 <- prepare(TE, seTE, treat1, treat2, studlab, func.inverse = func.inverse)
  ps <- prepare(TE, seTE, agent1, agent2, studlab, func.inverse = func.inverse)
  #
  o <- order(p0$order)
  os <- order(ps$order)
  
  #
  # Standard network meta-analysis
  #
  same_agents <- agent1 == agent2
  #
  Q.lump <- df.Q.lump <- pval.Q.lump <- NA
  #
  Q.split <- df.Q.split <- pval.Q.split <- NA
  #
  if (!any(same_agents)) {
    #
    # Lumping approach
    #
    if (netconnection(agent1, agent2, ps$studlab[os])$n.subnets == 1) {
      net1 <- netmeta(ps$TE[os], ps$seTE[os], agent1, agent2, ps$studlab[os],
                      reference.group = reference.group,
                      tol.multiarm = tol.multiarm,
                      tol.multiarm.se = tol.multiarm.se,
                      details.chkmultiarm = details.chkmultiarm
      )
      #
      Q.lump <- net1$Q
      df.Q.lump <- net1$df.Q
      pval.Q.lump <- net1$pval.Q
    }
    #
    # Splitting approach
    #
    if (netconnection(treat1, treat2, ps$studlab[os])$n.subnets == 1) {
      net2 <- netmeta(ps$TE[os], ps$seTE[os], treat1, treat2, ps$studlab[os],
                      reference.group = reference.group,
                      tol.multiarm = tol.multiarm,
                      tol.multiarm.se = tol.multiarm.se,
                      details.chkmultiarm = details.chkmultiarm
      )
      #
      Q.split <- net2$Q
      df.Q.split <- net2$df.Q
      pval.Q.split <- net2$pval.Q
    }
  }
  
  #
  # Common effects DR-NMA model
  #
  res.c <- nma_dose(p0$TE[o], p0$seTE[o], p0$weights[o], p0$studlab[o],
                    agent1, agent2, p0$treat1[o], p0$treat2[o],
                    p0$narms[o],
                    Xd = X.matrix, D = D.matrix, n = n,
                    level = level, reference = reference.group
  )
  #
  # Random effects DR-NMA model
  #
  tau <- res.c$tau
  #
  p1 <- prepare(TE, seTE, treat1, treat2, studlab, tau = tau, func.inverse = func.inverse)
  #
  res.r <- nma_dose(p1$TE[o], p1$seTE[o], p1$weights[o], p1$studlab[o],
                    agent1, agent2, p1$treat1[o], p1$treat2[o],
                    p1$narms[o],
                    Xd = X.matrix, D = D.matrix, n = n,
                    level = level, reference = reference.group)
  
  
  #
  #
  # (6) Generate DR-NMA object
  #
  #
  
  res <- list(
    studlab = studlab,
    agent1 = agent1,
    agent2 = agent2,
    dose1 = dose1,
    dose2 = dose2,
    treat1 = treat1,
    treat2 = treat2,
    #
    TE = p0$TE[o],
    seTE = p0$seTE[o],
    seTE.adj.common = sqrt(1 / p0$weights[o]),
    seTE.adj.random = sqrt(1 / p1$weights[o]),
    #
    k = res.c$k,
    m = res.c$m,
    a = res.c$a,
    n = n,
    trts = trts,
    agents = agents,
    inactive = inactive,
    #
    common.dose = common.dose,
    #
    TE.common = res.c$all.comparisons$TE,
    seTE.common = res.c$all.comparisons$seTE,
    lower.common = res.c$all.comparisons$lower,
    upper.common = res.c$all.comparisons$upper,
    statistic.common = res.c$all.comparisons$statistic,
    pval.common = res.c$all.comparisons$p,
    #
    TE.common = res.c$all.comparisons$TE,
    seTE.common = res.c$all.comparisons$seTE,
    lower.common = res.c$all.comparisons$lower,
    upper.common = res.c$all.comparisons$upper,
    statistic.common = res.c$all.comparisons$statistic,
    pval.common = res.c$all.comparisons$p,
    #
    TE.drnma.common = res.c$TE.nma,
    seTE.drnma.common = res.c$seTE.nma,
    lower.drnma.common = res.c$lower.nma,
    upper.drnma.common = res.c$upper.nma,
    statistic.drnma.common = res.c$statistic.nma,
    pval.drnma.common = res.c$pval.nma,
    #
    TE.random = res.r$all.comparisons$TE,
    seTE.random = res.r$all.comparisons$seTE,
    lower.random = res.r$all.comparisons$lower,
    upper.random = res.r$all.comparisons$upper,
    statistic.random = res.r$all.comparisons$statistic,
    pval.random = res.r$all.comparisons$p,
    #
    TE.drnma.random = res.r$TE.nma,
    seTE.drnma.random = res.r$seTE.nma,
    lower.drnma.random = res.r$lower.nma,
    upper.drnma.random = res.r$upper.nma,
    statistic.drnma.random = res.r$statistic.nma,
    pval.drnma.random = res.r$pval.nma,
    #
    Q = res.c$Q,
    df.Q = res.c$df.Q,
    pval.Q = res.c$pval.Q,
    Q.to.df.ratio = res.c$H^2,
    #
    tau = res.c$tau,
    tauml = res.c$tauML,
    I2 = res.c$I2, lower.I2 = res.c$lower.I2, upper.I2 = res.c$upper.I2,
    #
    Q.lump = Q.lump,
    df.Q.lump = df.Q.lump,
    pval.Q.lump = pval.Q.lump,
    #
    Q.split = Q.split,
    df.Q.split = df.Q.split,
    pval.Q.split = pval.Q.split,
    #
    B.matrix = B.matrix,
    D_obs.matrix = D_obs.matrix,
    D.matrix = D.matrix,
    X.matrix = X.matrix,
    #
    sm = sm,
    level = level,
    common = common,
    random = random,
    method = method,
    #
    param = param,
    knots = if (method == "rcs") attr(X.matrix, "param"),
    #
    reference.group = reference.group,
    #
    func.inverse = deparse(substitute(func.inverse)),
    #
    backtransf = backtransf,
    #
    data = if (keepdata) data else NULL
  )
  #
  class(res) <- "netdose"
  #
  res
}
