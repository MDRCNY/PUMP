# Helper function for Westfall Young

enforces monotonicity in p-values.

## Usage

``` r
get_adjp_minp(ind.B, rawp.order)
```

## Arguments

- ind.B:

  matrix of indicator variables for if each raw test statistic exceeds
  the null test statistics. dimensions: nrow = tnum, ncol = M.

- rawp.order:

  order of raw p-values in ascending order

## Value

returns adjusted p-value matrix
