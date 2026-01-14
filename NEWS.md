## netdose, version 0.7-4 (14-01-2025)

### Major changes

* netdose():
  - exponential model with a parameter p included

## netdose, version 0.7-2 

### Major changes

* Return Q statistic of standard network meta-analysis model using the splitting
  approach in addition to lumping approach

* Use R function facet_wrap() from the R package **ggplot2** for panels in
  dose-response plot

* New function dotplot() to produce a dot plot showing the dose distribution

* R packages **dplyr**, **tidyr** and **ggh4x** added to Imports

* R packages **grid**, **gridExtra** and **rlang** removed from Imports

### Bug fixes

* netdose():
  - remove test statistic for difference in goodness of fit between standard
    and dose-response network meta-analysis model

* predict.netdose():
  - calculate correct predictions for restricted cubic splines

* plot.netdose():
  - study results were not printed in some networks

### User-visible changes

* plot.netdose():
  - new defaults for arguments 'benchmark.threshold' and 'plateau.threshold'
  - new arguments 'col.line', 'col.bmdl' and 'col.med' to change the colour of
    dose-response curves and other information
  - new argument 'legend' to suppress printing of legend
  - argument 'same.ylim' removed (instead use same y-limits in panels if
    argument 'ylim' is not NULL)

### Internal changes

* netdose():
  - list elements 'Q.diff', 'df.Q.diff' and 'pval.Q.diff' removed
  - list elements 'Q.dose', 'df.Q.dose' and 'pval.Q.dose' renamed to
    'Q', 'df.Q' and 'pval.Q'
  - list elements 'Q.standard', 'df.Q.standard' and 'pval.Q.standard' renamed to
    'Q.lump', 'df.Q.lump' and 'pval.Q.lump'
  - new list elements 'Q.split', 'df.Q.split' and 'pval.Q.split'
  - rename list element 'Xd' to 'X.matrix'
  - new list element 'B.matrix'
  - new list element 'D_obs.matrix' with the observed doses
  - new list element 'D.matrix' with the transformed doses


## netdose, version 0.6-0

- initial release on GitHub

