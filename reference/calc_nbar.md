# This function calculates needed nbar to achieve a given power

This function calculates needed nbar to achieve a given power

## Usage

``` r
calc_nbar(
  d_m,
  MT = 2.8,
  MDES,
  J = NULL,
  K = NULL,
  Tbar,
  R2.1,
  R2.2 = NULL,
  ICC.2 = NULL,
  omega.2 = NULL,
  R2.3 = NULL,
  ICC.3 = NULL,
  omega.3 = NULL
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

- J:

  scalar; the number of schools

- K:

  scalar; the number of level 3 units (districts).

- Tbar:

  scalar; the proportion of samples that are assigned to the treatment.

- R2.1:

  scalar, or vector of length M; percent of variation explained by level
  1 covariates for each outcome.

- R2.2:

  scalar, or vector of length M; percent of variation explained by level
  2 covariates for each outcome.

- ICC.2:

  scalar, or vector of length M; level 2 (school) intraclass
  correlation.

- omega.2:

  scalar, or vector of length M; ratio of variance of level 2 average
  impacts to variance of level 2 random intercepts.

- R2.3:

  scalar, or vector of length M; percent of variation explained by level
  3 covariates for each outcome.

- ICC.3:

  scalar, or vector length M; level 3 (district) intraclass correlation.

- omega.3:

  scalar, or vector of length M; ratio of variance of level 3 average
  impacts to variance of level 3 random intercepts.

## Value

nbar, the number of individuals needed, or NA if not possible given
design
