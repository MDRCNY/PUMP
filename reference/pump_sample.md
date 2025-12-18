# Estimate the required sample size (core function)

The user chooses the context (d_m), MTP, type of sample size, MDES,
power definition, and choices of all relevant design parameters.

The functions performs a search algorithm, and returns the sample size
value within the specified tolerance. For a list of choices for specific
parameters, see pump_info().

## Usage

``` r
pump_sample(
  d_m,
  MTP = NULL,
  typesample,
  MDES,
  M = 1,
  numZero = NULL,
  nbar = NULL,
  J = NULL,
  K = NULL,
  target.power,
  power.definition,
  alpha,
  two.tailed = TRUE,
  Tbar,
  numCovar.1 = 0,
  numCovar.2 = 0,
  numCovar.3 = 0,
  R2.1 = 0,
  R2.2 = 0,
  R2.3 = 0,
  ICC.2 = 0,
  ICC.3 = 0,
  rho = NULL,
  rho.matrix = NULL,
  omega.2 = 0,
  omega.3 = 0,
  B = 1000,
  max.steps = 20,
  tnum = 1000,
  start.tnum = tnum/10,
  final.tnum = 4 * tnum,
  parallel.WY.cores = 1,
  updateProgress = NULL,
  max_sample_size_nbar = 10000,
  max_sample_size_JK = 1000,
  tol = 0.01,
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

- typesample:

  string; type of sample size to calculate: "nbar", "J", or "K".

- MDES:

  scalar or vector; the desired MDES values for each outcome. Please
  provide a scalar, a vector of length M, or vector of values for
  non-zero outcomes.

- M:

  scalar; the number of hypothesis tests (outcomes), including zero
  outcomes.

- numZero:

  scalar; additional number of outcomes assumed to be zero. Please
  provide NumZero + length(MDES) = M, if length(MDES) is not 1.

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

  target power for search algorithm.

- power.definition:

  see pump_info() for possible power definitions.

- alpha:

  scalar; the family wise error rate (FWER).

- two.tailed:

  scalar; TRUE/FALSE for two-tailed or one-tailed power calculation.

- Tbar:

  scalar; the proportion of samples that are assigned to the treatment.

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

- rho:

  scalar; assumed correlation between all pairs of test statistics.

- rho.matrix:

  matrix; alternate specification allowing a full matrix of correlations
  between test statistics. Must specify either rho or rho.matrix, but
  not both.

- omega.2:

  scalar, or vector of length M; ratio of variance of level 2 average
  impacts to variance of level 2 random intercepts.

- omega.3:

  scalar, or vector of length M; ratio of variance of level 3 average
  impacts to variance of level 3 random intercepts.

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

- max_sample_size_nbar:

  scalar; default upper bound for nbar for search algorithm.

- max_sample_size_JK:

  scalar; default upper bound for J or K for search algorithm.

- tol:

  tolerance for target power, defaults to 0.01 (1 This parameter
  controls when the search is done: when estimated power (checked with
  \`final.tnum\` iterations) is within \`tol\`, the search stops.

- give.optimizer.warnings:

  whether to return verbose optimizer warnings.

- verbose:

  TRUE/FALSE; Print out diagnostics of time, etc.

## Value

a pumpresult object containing sample size results.

## See also

For more detailed information about this function and the user choices,
see the manuscript \<doi:10.18637/jss.v108.i06\>, which includes a
detailed Technical Appendix including information about the designs and
models and parameters.

## Examples

``` r
J <- pump_sample(
  d_m = 'd2.1_m2fc',
  MTP = 'HO',
  power.definition = 'D1indiv',
  typesample = 'J',
  target.power = 0.8,
  nbar = 50,
  M = 3,
  MDES = 0.125,
  Tbar = 0.5, alpha = 0.05,
  numCovar.1 = 1,
  R2.1 = 0.1, ICC.2 = 0.05, rho = 0.2,
  tnum = 1000)
```
