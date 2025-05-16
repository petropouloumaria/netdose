#'  Antidepressants Dataset for Dose-Response Network Meta-Analysis
#'
#' @description
#' This dataset is a synthesis of randomized controlled studies investigating
#' the effects of antidepressant interventions for unipolar major depressive
#' disorder.
#'
#' @details
#' The dataset includes data on several antidepressant agents and their observed
#' effects across multiple clinical randomized controlled trials. The primary
#' outcome is the response rate, defined as the proportion of patients achieving
#' at least a 50% reduction in symptoms.
#'
#' The dataset is structured in an arm-level format and includes the following
#' variables:
#' \itemize{
#'   \item \code{drug}: Name of the antidepressant agent.
#'   \item \code{r}: Number of participants who responded to treatment.
#'   \item \code{n}: Total number of participants in the study arm.
#'   \item \code{dose}: Dose level of the antidepressant agent.
#'   \item \code{studyid}: Unique study identifier.
#' }
#'
#' This dataset is intended for use in dose-response network meta-analysis to
#' explore the effects of the several agents across various doses.
#'
#' @name antidepressants
#' @aliases antidepressants
#'
#' @docType data
#'
#' @format
#' A data frame with the following columns:
#' \tabular{rl}{
#'   \bold{\emph{drug}} \tab Character vector indicating the name of the
#'     antidepressant agent. \cr
#'   \bold{\emph{r}} \tab Integer vector for the number of participants who
#'     responded to treatment. \cr
#'   \bold{\emph{n}} \tab Integer vector for the total number of participants
#'     in the study arm. \cr
#'   \bold{\emph{dose}} \tab Numeric vector specifying the dose level of the
#'     antidepressant agent. \cr
#'   \bold{\emph{studyid}} \tab Character vector with unique study identifiers.
#' }
#'
#' @source
#' The dataset is a subset of data derived from:
#' 
#' Cipriani A et al. (2018):
#' Comparative efficacy and acceptability of 21 antidepressant drugs for the
#' acute treatment of adults with major depressive disorder: a systematic review
#' and network meta-analysis.
#' \emph{Lancet},
#' \bold{391}, 1357--66
#' 
#' @keywords datasets
#'
#' @examples
#' summary(antidepressants)

"antidepressants"
