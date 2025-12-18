# Calculating Needed Sample Size for Raw (Unadjusted) Power

This is a helper function for getting a needed Sample Size when no
adjustments has been made to the test statistics.

## Usage

``` r
pump_sample_raw(
  d_m,
  MTP,
  typesample,
  MDES,
  nbar = NULL,
  J = NULL,
  K = NULL,
  target.power,
  Tbar,
  alpha = 0.05,
  two.tailed,
  numCovar.1 = 0,
  numCovar.2 = 0,
  numCovar.3 = 0,
  R2.1,
  R2.2 = NULL,
  R2.3 = NULL,
  ICC.2 = NULL,
  ICC.3 = NULL,
  omega.2 = NULL,
  omega.3 = NULL,
  max.steps = 100,
  warn.small = FALSE
)
```

## Arguments

- d_m:

  string; a single context, which is a design and model code. See
  pump_info() for list of choices.

- MTP:

  string, or vector of strings; multiple testing procedure(s). See
  pump_info() for list of choices.

- typesample:

  type of sample size to calculate: J, K, or nbar

- MDES:

  scalar or vector; the desired MDES values for each outcome. Please
  provide a scalar, a vector of length M, or vector of values for
  non-zero outcomes.

- nbar:

  scalar; the harmonic mean of the number of level 1 units per level 2
  unit (students per school). Note that this is not the total number of
  level 1 units, but instead the number of level 1 units nested within
  each level 2 unit, so the total number of level 1 units is nbar x J x
  K.

- J:

  scalar; the harmonic mean of number of level 2 units per level 3 unit
  (schools per district). Note that this is not the total number of
  level 2 units, but instead the number of level 2 units nested within
  each level 3 unit, so the total number of level 2 units is J x K.

- K:

  scalar; the number of level 3 units (districts).

- target.power:

  target power to arrive at

- Tbar:

  scalar; the proportion of samples that are assigned to the treatment.

- alpha:

  scalar; the family wise error rate (FWER).

- two.tailed:

  scalar; TRUE/FALSE for two-tailed or one-tailed power calculation.

- numCovar.1:

  scalar; number of level 1 (individual) covariates.

- numCovar.2:

  scalar; number of level 2 (school) covariates.

- numCovar.3:

  scalar; number of level 3 (district) covariates.

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

- max.steps:

  how many steps allowed before terminating

- warn.small:

  Warn if degrees of freedom issues are causing inability to achieve
  target power for sample size.

## Value

Requisite sample size (as integer) and associated degrees of freedom.

## Details

It is for a single, individual outcome. It only takes scalar values for
all its arguments, and does not have an M argument (for number of
outcomes).

It requires iteration because we do not know the degrees of freedom, and
so we guess and then calculate sample size, and then recalculate df
based on sample size, until we converge.

It is possible that the returned sample size will be the minimum sample
size required to have at least 1 degree of freedom (even if this
provides higher than target level power).
