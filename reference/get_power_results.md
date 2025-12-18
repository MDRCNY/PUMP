# Calculates different definitions of power (support function)

This function takes in a matrix of adjusted p-values and unadjusted
p-values and outputs different types of power.

This function is mostly for internal use, but may be of interest to
users who wish to calculate power on their own.

## Usage

``` r
get_power_results(
  adj.pval.mat,
  unadj.pval.mat,
  ind.nonzero,
  alpha,
  drop.zero.outcomes = TRUE,
  adj = TRUE
)
```

## Arguments

- adj.pval.mat:

  matrix; adjusted p-values, columns are outcomes

- unadj.pval.mat:

  matrix; unadjusted p-values, columns are outcomes

- ind.nonzero:

  vector; which outcomes correspond to nonzero effects.

- alpha:

  scalar; the family wise error rate (FWER).

- drop.zero.outcomes:

  logical; whether to report power results for outcomes with MDES = 0.

- adj:

  logical; whether p-values are unadjusted or not.

## Value

data frame; power results for individual, minimum, complete power.
