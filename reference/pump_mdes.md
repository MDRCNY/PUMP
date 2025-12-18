# Estimate the minimum detectable effect size (MDES) (core function)

The user chooses the context (d_m), MTP, power definition, and choices
of all relevant design parameters.

The functions performs a search algorithm, and returns the MDES value
within the specified tolerance. For a list of choices for specific
parameters, see pump_info().

## Usage

``` r
pump_mdes(
  d_m,
  MTP = NULL,
  numZero = NULL,
  propZero = NULL,
  M = 1,
  nbar,
  J = 1,
  K = 1,
  Tbar,
  alpha = 0.05,
  two.tailed = TRUE,
  target.power = 0.8,
  power.definition,
  tol = 0.02,
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
  rho = NULL,
  rho.matrix = NULL,
  B = 1000,
  max.steps = 20,
  tnum = 1000,
  start.tnum = round(tnum/10),
  final.tnum = 4 * tnum,
  parallel.WY.cores = 1,
  updateProgress = NULL,
  give.optimizer.warnings = FALSE,
  verbose = FALSE
)
```

## Arguments

- d_m:

  string; a single context, which is a design and model code. See
  pump_info() for list of choices.

- MTP:

  string, or vector of strings; multiple testing procedure(s). See
  pump_info() for list of choices.

- numZero:

  scalar; additional number of outcomes assumed to be zero. Please
  provide NumZero + length(MDES) = M, if length(MDES) is not 1.

- propZero:

  scalar; proportion of outcomes assumed to be zero (alternative
  specification to numZero). length(MDES) should be 1 or equal to
  (1-propZero)\*M.

- M:

  scalar; the number of hypothesis tests (outcomes), including zero
  outcomes.

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

- Tbar:

  scalar; the proportion of samples that are assigned to the treatment.

- alpha:

  scalar; the family wise error rate (FWER).

- two.tailed:

  scalar; TRUE/FALSE for two-tailed or one-tailed power calculation.

- target.power:

  target power for search algorithm.

- power.definition:

  see pump_info() for possible power definitions.

- tol:

  tolerance for target power, defaults to 0.01 (1 This parameter
  controls when the search is done: when estimated power (checked with
  \`final.tnum\` iterations) is within \`tol\`, the search stops.

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

- rho.matrix:

  matrix; alternate specification allowing a full matrix of correlations
  between test statistics. Must specify either rho or rho.matrix, but
  not both.

- B:

  scalar; the number of permutations for Westfall-Young procedures.

- max.steps:

  how many steps allowed before terminating.

- tnum:

  max number of samples for first iteration of search algorithm.

- start.tnum:

  number of samples to start search (this will increase with each step).

- final.tnum:

  number of samples for final draw.

- parallel.WY.cores:

  number of cores to use for parallel processing of WY-SD.

- updateProgress:

  function to update progress bar (only used for PUMP shiny app).

- give.optimizer.warnings:

  whether to return verbose optimizer warnings.

- verbose:

  TRUE/FALSE; Print out diagnostics of time, etc.

## Value

a pumpresult object containing MDES results.

## See also

For more detailed information about this function and the user choices,
see the manuscript \<doi:10.18637/jss.v108.i06\>, which includes a
detailed Technical Appendix including information about the designs and
models and parameters.

## Examples

``` r
mdes <-  pump_mdes(
  d_m = "d3.1_m3rr2rr",
  MTP = 'HO',
  power.definition = 'D1indiv',
  target.power = 0.6,
  J = 30,
  K = 15,
  nbar = 50,
  M = 3,
  Tbar = 0.5, alpha = 0.05,
  two.tailed = FALSE,
  numCovar.1 = 1, numCovar.2 = 1,
  R2.1 = 0.1, R2.2 = 0.1,
  ICC.2 = 0.2, ICC.3 = 0.2,
  omega.2 = 0.1, omega.3 = 0.1,
  rho = 0.5, tnum = 2000)
```
