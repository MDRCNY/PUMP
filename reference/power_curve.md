# Obtain a power curve for a range of sample size or MDES values

This is used to see how power changes as a function of sample size or
MDES. It takes a fit pumpresult and calculates a power curve based on
that scenario coupled with a passed range of values to make the curve
over.

## Usage

``` r
power_curve(
  x,
  all = FALSE,
  low = NULL,
  high = NULL,
  grid.size = 5,
  tnum = 2000
)
```

## Arguments

- x:

  a pumpresult object.

- all:

  logical; if TRUE, merge in the search path from the original search.

- low:

  scalar; low range for curve.

- high:

  scalar; high range for the curve.

- grid.size:

  scalar; number of points to calculate power for.

- tnum:

  scalar; number of iterations to calculate power at each grid point.

## Value

data.frame of power results.
