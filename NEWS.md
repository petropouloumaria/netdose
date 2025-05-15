## netdose, version 0.7-0 (2025-mm-dd)

### Major changes

* Return Q statistic of standard network meta-analysis model using the splitting
  approach in addition to lumping approach

### Bug fixes

* netdose():
  - remove test statistic for difference in goodness of fit between standard
    and dose-response model

* plot.netdose():
  - no study results printed in some networks

### Internal changes

* netdose():
  - list elements 'Q.diff', 'df.Q.diff' and 'pval.Q.diff' removed
  - list elements 'Q.dose', 'df.Q.dose' and 'pval.Q.dose' renamed to
    'Q', 'df.Q' and 'pval.Q'
  - list elements 'Q.standard', 'df.Q.standard' and 'pval.Q.standard' renamed to
    'Q.lump', 'df.Q.lump' and 'pval.Q.lump'
  - new list elements 'Q.split', 'df.Q.split' and 'pval.Q.split'


## netdose, version 0.6-0

- initial release on GitHub

