netdose: Dose-Response Network Meta-Analysis in a Frequentist Way

[![License: GPL
(\>=2)](https://img.shields.io/badge/license-GPL-blue)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)
[![CRAN
Version](https://www.r-pkg.org/badges/version/netdose)](https://cran.r-project.org/package=netdose)
[![CRAN
Downloads](https://cranlogs.r-pkg.org/badges/netdose)](https://cranlogs.r-pkg.org/badges/netdose)

## Description

`netdose` is an R package for conducting **Dose-Response Network
Meta-Analysis (DR-NMA)** in the frequentist framework. It supports
multiple dose-response functions, including linear, exponential,
fractional polynomials, and restricted cubic splines.

## Installation

You can install the **netdose** package from GitHub repository as
follows:

Installation using R package
**[remotes](https://cran.r-project.org/package=remotes)**:

```r
install.packages("remotes")
remotes::install_github("petropouloumaria/netdose")
```

## Usage

You can load the **netdose** package.

```r
library("netdose")
```

Use the anesthesia dataset from **netdose** package. This dataset is a
synthesis of randomized controlled studies investigating the effects of
different agents for preventing vomiting in adults after general
anaesthesia. The outcome is the occurrence of vomiting within 24 hours
after surgery. Multicomponent interventions have been excluded from this
subset to focus on single-agent interventions and their dose-response
relationships.

References:

Weibel S, Rücker G, Eberhart LHJ, et al. *Drugs for preventing
postoperative nausea and vomiting in adults after general anaesthesia: a
network meta-analysis.*Cochrane Database of Systematic Reviews.
2020;**10**(10):CD012859.

Weibel S, Schaefer MS, Raj D, et al. *Drugs for preventing postoperative
nausea and vomiting in adults after general anaesthesia: an abridged
Cochrane network meta-analysis.* Anaesthesia. 2021;**76**:962-973.

Transform data from arm-based to contrast-based format using the
function **pairwise** from **meta** package.

```r
dat <- pairwise(
  agent = list(agent1, agent2, agent3, agent4, agent5),
  event = list(event1, event2, event3, event4, event5),
  n = list(n1, n2, n3, n4, n5),
  dose = list(dose1, dose2, dose3, dose4, dose5),
  data = anesthesia,
  studlab = study
)
```

Perform the DR-NMA model with different dose-response functions.

Linear Dose-Response NMA Model

```r
dr_l <- netdose(TE, seTE, agent1, dose1, agent2, dose2, 
   studlab, data = dat)
```

Exponential Dose-Response NMA Model

```r
dr_exp <- netdose(TE, seTE, agent1, dose1, agent2, dose2, 
   studlab, data = dat, method = "exponential")
```

Fractional Polynomial (FP1) with power p = 0.5

```r
dr_fp1_p05 <- netdose(TE, seTE, agent1, dose1, agent2, dose2, 
   studlab, data = dat, method = "fp1", param = 0.5)
```

Restricted Cubic Splines (RCS) with knots k = (0.1, 0.5, 0.9)

```r
dr_rcs <- netdose(TE, seTE, agent1, dose1, agent2, dose2, 
   studlab, data = dat, method = "rcs")
```

Predicting Dose-Response Relationships

Predicted Values for Different DR-NMA Models (exponential, FP1 (p =
0.5), RCS (k = (0.1, 0.5, 0.9)))

```r
predict(dr_exp)
predict(dr_fp1_p05)
predict(dr_rcs)
```

Visualizing Dose-Response Relationships

``` r
plot(dr_exp, only.direct = FALSE)
plot(dr_fp1_p05, only.direct = FALSE)
plot(dr_rcs, only.direct = FALSE)
```

## License

This package is released under the **GPL-2 License**.
