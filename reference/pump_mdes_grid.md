# Run pump_mdes on varying values of parameters (grid function)

See pump_power_grid() for more details.

## Usage

``` r
pump_mdes_grid(
  d_m,
  MTP = NULL,
  M = 1,
  target.power = 0.8,
  power.definition = NULL,
  tol = 0.01,
  propZero = NULL,
  numZero = NULL,
  nbar,
  J = 1,
  K = 1,
  Tbar = 0.5,
  alpha = 0.05,
  numCovar.1 = NULL,
  numCovar.2 = NULL,
  numCovar.3 = NULL,
  R2.1 = NULL,
  R2.2 = NULL,
  R2.3 = NULL,
  ICC.2 = NULL,
  ICC.3 = NULL,
  omega.2 = NULL,
  omega.3 = NULL,
  rho = NULL,
  verbose = FALSE,
  drop.unique.columns = TRUE,
  ...
)
```

## Arguments

- d_m:

  string; a single context, which is a design and model code. See
  pump_info() for list of choices.

- MTP:

  string, or vector of strings; multiple testing procedure(s). See
  pump_info() for list of choices.

- M:

  scalar; the number of hypothesis tests (outcomes), including zero
  outcomes.

- target.power:

  target power for search algorithm.

- power.definition:

  see pump_info() for possible power definitions.

- tol:

  tolerance for target power, defaults to 0.01 (1 This parameter
  controls when the search is done: when estimated power (checked with
  \`final.tnum\` iterations) is within \`tol\`, the search stops.

- propZero:

  scalar; proportion of outcomes assumed to be zero (alternative
  specification to numZero). length(MDES) should be 1 or equal to
  (1-propZero)\*M.

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

- Tbar:

  scalar; the proportion of samples that are assigned to the treatment.

- alpha:

  scalar; the family wise error rate (FWER).

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

- verbose:

  TRUE/FALSE; Print out diagnostics of time, etc.

- drop.unique.columns:

  logical; drop all parameter columns that did not vary across the grid.

- ...:

  extra arguments passed to the underlying pump_power, pump_sample, or
  pump_mdes functions.

## Value

a pumpgridresult object containing MDES results.

## See also

Other grid functions:
[`pump_power_grid()`](https://mdrcny.github.io/PUMP/reference/pump_power_grid.md),
[`pump_sample_grid()`](https://mdrcny.github.io/PUMP/reference/pump_sample_grid.md)

## Examples

``` r
g <- pump_mdes_grid(d_m = "d3.2_m3ff2rc", MTP = "HO",
  target.power = c( 0.50, 0.80 ), power.definition = "D1indiv",
  tol = 0.05, M = 5, J = c( 3, 9 ), K = 7, nbar = 58,
  Tbar = 0.50, alpha = 0.15, numCovar.1 = 1, numCovar.2 = 1,
  R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, ICC.3 = 0.9,
  rho = 0.4, tnum = 200)
```
