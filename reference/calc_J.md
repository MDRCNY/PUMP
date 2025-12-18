# This function calculates needed J to achieve a given (unadjusted) power

This function calculates needed J to achieve a given (unadjusted) power

## Usage

``` r
calc_J(
  d_m,
  MT = 2.8,
  MDES,
  K = NULL,
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

  a single RCT design (see list/naming convention)

- MT:

  Number of approximate effect-size unit SEs (adjusted for degrees of
  freedom issues) that the MDES needs to be to achieve desired power.
  E.g., 2.8 for normal theory.

- MDES:

  scalar; the MDES values for each outcome

- K:

  scalar; the number of level 3 units (districts).

- nbar:

  scalar; the harmonic mean of the number of level 1 units per level 2
  unit (students per school). Note that this is not the total number of
  level 1 units, but instead the number of level 1 units nested within
  each level 2 unit, so the total number of level 1 units is nbar x J x
  K.

- Tbar:

  scalar; the proportion of samples that are assigned to the treatment.

- R2.1:

  scalar, or vector of length M; percent of variation explained by level
  1 covariates for each outcome.

- R2.2:

  scalar, or vector of length M; percent of variation explained by level
  2 covariates for each outcome.

- R2.3:

  scalar, or vector of length M; percent of variation explained by level
  3 covariates for each outcome.

- ICC.2:

  scalar, or vector of length M; level 2 (school) intraclass
  correlation.

- ICC.3:

  scalar, or vector length M; level 3 (district) intraclass correlation.

- omega.2:

  scalar, or vector of length M; ratio of variance of level 2 average
  impacts to variance of level 2 random intercepts.

- omega.3:

  scalar, or vector of length M; ratio of variance of level 3 average
  impacts to variance of level 3 random intercepts.

## Value

J, the number of schools needed
