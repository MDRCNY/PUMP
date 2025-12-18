# Computes Q_m, the standard error of the effect size estimate

Function to calculate the theoretical true (unadjusted) standard error
of the ATE estimate for a given d_m and model, in effect size units.

## Usage

``` r
calc_SE(
  d_m,
  J,
  K,
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

  a single RCT d_m (see list/naming convention).

- J:

  scalar; the number of schools

- K:

  scalar; the number of districts

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

  scalar, or vector of length M; school intraclass correlation

- ICC.3:

  scalar, or vector of length M; district intraclass correlation

- omega.2:

  scalar, or vector of length M; ratio of school effect size variability
  to random effects variability

- omega.3:

  scalar, or vector of length M; ratio of district effect size
  variability to random effects variability

## Value

vector; the standard error of the effect size estimate
