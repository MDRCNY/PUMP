# Westfall Young Step Down Function

This adjustment function utilizes the comp_rawp_ss helper function to
compare each row of the matrix sample p-values under alternative
hypothesis to all the rows in the matrix of the p-values under the
complete null.

## Usage

``` r
adjp_wysd(
  rawp.mat,
  B,
  Sigma,
  t.df,
  two.tailed,
  cl = NULL,
  verbose = TRUE,
  updateProgress = NULL
)
```

## Arguments

- rawp.mat:

  a matrix of raw p-values under H1. dimension: nrow = tnum, ncol = M

- B:

  numer of WY permutations

- Sigma:

  correlation matrix of null p-values

- t.df:

  degrees of freedom of null p-values

- two.tailed:

  one or two-tailed test

- cl:

  cluster object for parallel computing

- verbose:

  whether to print out messaging

- updateProgress:

  function to update progress bar (only used for PUMP shiny app)

## Value

a matrix of adjusted p-values
