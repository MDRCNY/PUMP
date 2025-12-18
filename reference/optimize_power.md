# Optimizes power to help in search for MDES or SS

Optimizes power to help in search for MDES or SS

## Usage

``` r
optimize_power(
  d_m,
  search.type,
  MTP,
  target.power,
  power.definition,
  tol,
  start.low,
  start.high,
  MDES = NULL,
  J = NULL,
  K = 1,
  nbar = NULL,
  M = M,
  numZero = numZero,
  Tbar = Tbar,
  alpha,
  two.tailed,
  numCovar.1 = 0,
  numCovar.2 = 0,
  numCovar.3 = 0,
  R2.1 = 0,
  R2.2 = 0,
  R2.3 = 0,
  ICC.2 = 0,
  ICC.3 = 0,
  omega.2 = 0,
  omega.3 = 0,
  rho,
  B = NULL,
  parallel.WY.cores = 1,
  max.steps = 20,
  tnum = 1000,
  start.tnum = round(tnum/10),
  final.tnum = 4 * tnum,
  give.warnings = FALSE,
  grid.only = FALSE,
  grid.size = 5
)
```

## Arguments

- d_m:

  string; a single context, which is a design and model code. See
  pump_info() for list of choices.

- search.type:

  type of optimization search being conducted. Options: K, J, nbar, mdes

- MTP:

  string, or vector of strings; multiple testing procedure(s). See
  pump_info() for list of choices.

- target.power:

  target power for search algorithm.

- power.definition:

  see pump_info() for possible power definitions.

- tol:

  tolerance for target power, defaults to 0.01 (1 This parameter
  controls when the search is done: when estimated power (checked with
  \`final.tnum\` iterations) is within \`tol\`, the search stops.

- start.low:

  lower bound for optimization procedure

- start.high:

  upper bound for optimization procedure

- MDES:

  scalar or vector; the desired MDES values for each outcome. Please
  provide a scalar, a vector of length M, or vector of values for
  non-zero outcomes.

- J:

  scalar; the harmonic mean of number of level 2 units per level 3 unit
  (schools per district). Note that this is not the total number of
  level 2 units, but instead the number of level 2 units nested within
  each level 3 unit, so the total number of level 2 units is J x K.

- K:

  scalar; the number of level 3 units (districts).

- nbar:

  scalar; the harmonic mean of the number of level 1 units per level 2
  unit (students per school). Note that this is not the total number of
  level 1 units, but instead the number of level 1 units nested within
  each level 2 unit, so the total number of level 1 units is nbar x J x
  K.

- M:

  scalar; the number of hypothesis tests (outcomes), including zero
  outcomes.

- numZero:

  scalar; additional number of outcomes assumed to be zero. Please
  provide NumZero + length(MDES) = M, if length(MDES) is not 1.

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

- rho:

  scalar; assumed correlation between all pairs of test statistics.

- B:

  scalar; the number of permutations for Westfall-Young procedures.

- parallel.WY.cores:

  number of cores to use for parallel processing of WY-SD.

- max.steps:

  how many steps allowed before terminating.

- tnum:

  scalar; the number of test statistics to draw. Increasing tnum
  increases precision and computation time.

- start.tnum:

  number of samples to start search (this will increase with each step).

- final.tnum:

  number of samples for final draw.

- give.warnings:

  whether to return optimizer warnings

- grid.only:

  TRUE means generate a grid from start.low to start.high, but do not do
  iterative search. (Useful for mapping out the power curve rather than
  identifying a point of particular power.)

- grid.size:

  Number of points to check in initial search grid. Grid will be spaced
  as a quadratic sequence (e.g., 0, 1, 4, 9, 16 for a 0-16 span).

## Value

power
