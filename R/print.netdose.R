#' @method print netdose
#' @export


print.netdose <- function(x,
                          common = x$common,
                          random = x$random,
                          backtransf = x$backtransf,
                          ##
                          digits = gs("digits"),
                          digits.stat = gs("digits.stat"),
                          digits.pval = gs("digits.pval"),
                          digits.pval.Q = max(gs("digits.pval.Q"), 2),
                          digits.Q = gs("digits.Q"),
                          digits.tau2 = gs("digits.tau2"),
                          digits.tau = gs("digits.tau"),
                          digits.I2 = gs("digits.I2"),
                          ##
                          scientific.pval = gs("scientific.pval"),
                          ##
                          big.mark = gs("big.mark"),
                          ##
                          text.tau2 = gs("text.tau2"),
                          text.tau = gs("text.tau"),
                          text.I2 = gs("text.I2"),
                          ...) {
    ## Check class
    ##
    chkclass(x, "netdose")


    cat("Number of studies: k = ", x$k, "\n", sep = "")
    cat("Number of pairwise comparisons: m = ", x$m, "\n", sep = "")
    cat("Number of treatments: n = ", x$n, "\n", sep = "")
    cat("Number of agents: a = ", x$a, "\n", sep = "")

    ##
    cat("\n")

    ##
    ## (a) Results for comparisons
    ##
    sm <- x$sm
    reference.group <- x$reference.group
    baseline.reference <- x$baseline.reference

    if (!backtransf & (is_relative_effect(sm) | sm == "VE")) {
        sm.lab <- paste0("log", if (sm == "VE") "VR" else sm)
    } else {
        sm.lab <- sm
    }
    ##
    ci.lab <- paste0(round(100 * x$level, 1), "%-CI")
    ##
    TE.common <- x$TE.drnma.common
    lowTE.common <- x$lower.drnma.common
    uppTE.common <- x$upper.drnma.common
    statistic.common <- x$statistic.drnma.common
    pval.common <- x$pval.drnma.common
    ##
    TE.random <- x$TE.drnma.random
    lowTE.random <- x$lower.drnma.random
    uppTE.random <- x$upper.drnma.random
    statistic.random <- x$statistic.drnma.random
    pval.random <- x$pval.drnma.random
    ##
    comptext <- paste("comparison: agent vs", reference.group)

    ##
    ##
    ##
    noeffect <- 1L * (backtransf & is_relative_effect(sm))
    #
    #
    if (backtransf) {
        TE.common <- backtransf(TE.common, sm)
        lowTE.common <- backtransf(lowTE.common, sm)
        uppTE.common <- backtransf(uppTE.common, sm)
        ##
        TE.random <- backtransf(TE.random, sm)
        lowTE.random <- backtransf(lowTE.random, sm)
        uppTE.random <- backtransf(uppTE.random, sm)
        #
        # Switch lower and upper limit for VE if results have been
        # backtransformed
        #
        if (sm == "VE") {
            tmp.l <- lowTE.common
            lowTE.common <- uppTE.common
            uppTE.common <- tmp.l
            #
            tmp.l <- lowTE.random
            lowTE.random <- uppTE.random
            uppTE.random <- tmp.l
        }
    }
    ##
    TE.common <- round(TE.common, digits)
    lowTE.common <- round(lowTE.common, digits)
    uppTE.common <- round(uppTE.common, digits)
    statistic.common <- round(statistic.common, digits.stat)
    ##
    TE.random <- round(TE.random, digits)
    lowTE.random <- round(lowTE.random, digits)
    uppTE.random <- round(uppTE.random, digits)
    statistic.random <- round(statistic.random, digits.stat)
    ##
    dat1.c <-
        cbind(
            formatN(TE.common, digits,
                text.NA = "NA",
                big.mark = big.mark
            ),
            formatCI(
                formatN(round(lowTE.common, digits),
                    digits, "NA",
                    big.mark = big.mark
                ),
                formatN(round(uppTE.common, digits),
                    digits, "NA",
                    big.mark = big.mark
                )
            ),
            formatN(statistic.common, digits.stat,
                text.NA = "NA",
                big.mark = big.mark
            ),
            formatPT(pval.common,
                digits = digits.pval,
                scientific = scientific.pval
            )
        )

    trts <- rownames(dat1.c)

    dat1.c[trts == reference.group, ] <- rep(".", ncol(dat1.c))
    dimnames(dat1.c) <-
        list(trts, c(sm.lab, ci.lab, "z", "p-value"))
    ##
    ## if (!is.na(x$TE.drnma.common[trts == reference.group]) &&
    ##    x$TE.drnma.common[trts == reference.group] == noeffect) {
    ##    dat1.c[trts == reference.group, ] <- "."
    ## }
    # rownames(dat1.c) <- trts.abbr



    dat1.r <-
        cbind(
            formatN(TE.random, digits,
                text.NA = "NA",
                big.mark = big.mark
            ),
            formatCI(
                formatN(round(lowTE.random, digits),
                    digits, "NA",
                    big.mark = big.mark
                ),
                formatN(round(uppTE.random, digits),
                    digits, "NA",
                    big.mark = big.mark
                )
            ),
            formatN(statistic.random, digits.stat,
                text.NA = "NA",
                big.mark = big.mark
            ),
            formatPT(pval.random,
                digits = digits.pval,
                scientific = scientific.pval
            )
        )

    dat1.r[trts == reference.group, ] <- rep(".", ncol(dat1.r))
    dimnames(dat1.r) <-
        list(trts, c(sm.lab, ci.lab, "z", "p-value"))
    ##
    ## if (!is.na(x$TE.drnma.random[trts == reference.group]) &&
    ##    x$TE.drnma.random[trts == reference.group] == noeffect) {
    ##    dat1.r[trts == reference.group, ] <- "."
    ## }
    # rownames(dat1.c) <- trts.abbr



    if (common) {
        if (reference.group != "") {
            cat("Common effects model")
            cat("\n")
            ##
            cat("Treatment estimate (sm = '", sm.lab,
                "', ", comptext, "):\n",
                sep = ""
            )
            prmatrix(dat1.c, quote = FALSE, right = TRUE)
            cat("\n")
        }
    }
    if (random) {
        if (reference.group != "") {
            cat("Random effects model")
            cat("\n")
            ##
            cat("Treatment estimate (sm = '", sm.lab,
                "', ", comptext, "):\n",
                sep = ""
            )
            prmatrix(dat1.r, quote = FALSE, right = TRUE)
            cat("\n")
        }
    }
    ##
    ## (d) Heterogeneity / inconsistency
    ##
    cat("Quantifying heterogeneity / inconsistency:\n",
        formatPT(x$tau^2,
            lab = TRUE, labval = text.tau2,
            digits = digits.tau2,
            lab.NA = "NA", big.mark = big.mark
        ),
        "; ",
        formatPT(x$tau,
            lab = TRUE, labval = text.tau,
            digits = digits.tau,
            lab.NA = "NA", big.mark = big.mark
        ),
        if (!is.na(x$I2)) {
            paste0("; ", text.I2, " = ", round(x$I2 * 100, digits.I2), "%")
        },
        if (!(is.na(x$lower.I2) | is.na(x$upper.I2))) {
            pasteCI(x$lower.I2 * 100, x$upper.I2 * 100, digits.I2, big.mark, unit = "%")
        },
        "\n\n",
        sep = ""
    )

    cat("Heterogeneity statistics:\n")
    ##
    hetdat <-
        data.frame(
            Q = formatN(
                c(
                    x$Q.dose,
                    x$Q.standard,
                    x$Q.diff
                ),
                digits.Q
            ),
            df.Q = formatN(c(
                x$df.Q.dose,
                x$df.Q.standard,
                x$df.Q.diff
            ), 0),
            pval = formatPT(
                c(
                    x$pval.Q.dose,
                    x$pval.Q.standard,
                    x$pval.Q.diff
                ),
                digits = digits.pval.Q,
                scientific = scientific.pval
            ),
            row.names = c(
                "DR-NMA model", "Standard model",
                "Difference"
            )
        )
    ##
    names(hetdat) <- c("Q", "df", "p-value")
    ##
    print(hetdat)

    invisible(NULL)
}
