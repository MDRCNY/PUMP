# Estimate power across definitions (core function)

The user chooses the context (d_m), MTP, MDES, and choices of all
relevant design parameters.

The functions returns power for all definitions of power for any MTP.
For a list of choices for specific parameters, see pump_info().

## Usage

``` r
pump_power(
  d_m,
  MTP = NULL,
  MDES,
  numZero = NULL,
  propZero = NULL,
  M = 1,
  nbar,
  J = 1,
  K = 1,
  Tbar,
  alpha = 0.05,
  two.tailed = TRUE,
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
  tnum = 10000,
  B = 1000,
  parallel.WY.cores = 1,
  drop.zero.outcomes = TRUE,
  updateProgress = NULL,
  validate.inputs = TRUE,
  long.table = FALSE,
  verbose = FALSE,
  exact.where.possible = TRUE
)
```

## Arguments

- d_m:

  string; a single context, which is a design and model code. See
  pump_info() for list of choices.

- MTP:

  string, or vector of strings; multiple testing procedure(s). See
  pump_info() for list of choices.

- MDES:

  scalar or vector; the desired MDES values for each outcome. Please
  provide a scalar, a vector of length M, or vector of values for
  non-zero outcomes.

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

- tnum:

  scalar; the number of test statistics to draw. Increasing tnum
  increases precision and computation time.

- B:

  scalar; the number of permutations for Westfall-Young procedures.

- parallel.WY.cores:

  number of cores to use for parallel processing of WY-SD.

- drop.zero.outcomes:

  whether to report power results for outcomes with MDES = 0. If ALL
  MDES = 0, then the first outcome will not be dropped.

- updateProgress:

  function to update progress bar (only used for PUMP shiny app).

- validate.inputs:

  TRUE/FALSE; whether or not to check whether parameters are valid given
  the choice of d_m.

- long.table:

  TRUE for table with power as rows, correction as columns, and with
  more verbose names. See \`transpose_power_table\`.

- verbose:

  TRUE/FALSE; Print out diagnostics of time, etc.

- exact.where.possible:

  TRUE/FALSE; whether to do exact calculations when M=1, or use
  simulation. Default is TRUE.

## Value

a pumpresult object containing power results.

## See also

For more detailed information about this function and the user choices,
see the manuscript \<doi:10.18637/jss.v108.i06\>, which includes a
detailed Technical Appendix including information about the designs and
models and parameters.

## Examples

``` r
pp <- pump_power(
   d_m = "d3.2_m3ff2rc",
   MTP = 'HO',
   nbar = 50,
   J = 30,
   K = 10,
   M = 5,
   MDES = 0.125,
   Tbar = 0.5, alpha = 0.05,
   numCovar.1 = 1, numCovar.2 = 1,
   R2.1 = 0.1, R2.2 = 0.1,
   ICC.2 = 0.2, ICC.3 = 0.2,
   omega.2 = 0, omega.3 = 0.1,
   rho = 0.5)
```
