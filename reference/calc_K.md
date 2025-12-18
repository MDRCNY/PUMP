# Calculates K, the number of districts

Calculates K, the number of districts

## Usage

``` r
calc_K(
  d_m,
  MT = 2.8,
  MDES,
  J,
  nbar,
  Tbar,
  R2.1,
  R2.2,
  R2.3,
  ICC.2,
  ICC.3,
  omega.2,
  omega.3
)
```

## Arguments

- d_m:

  a single RCT d_m (see list/naming convention)

- MT:

  multiplier

- MDES:

  scalar; the MDES value for all outcomes

- J:

  scalar; the number of schools

- nbar:

  scalar; the harmonic mean of the number of units per school

- Tbar:

  scalar; the proportion of samples that are assigned to the treatment

- R2.1:

  scalar, or vector of length M; percent of variation explained by Level
  1 covariates for each outcome

- R2.2:

  scalar, or vector of length M; percent of variation explained by Level
  2 covariates for each outcome

- R2.3:

  scalar, or vector of length M; percent of variation explained by Level
  3 covariates for each outcome

- ICC.2:

  scalar; school intraclass correlation

- ICC.3:

  scalar; district intraclass correlation

- omega.2:

  scalar; ratio of school effect size variability to random effects
  variability

- omega.3:

  scalar; ratio of district effect size variability to random effects
  variability

## Value

K, the number of districts
