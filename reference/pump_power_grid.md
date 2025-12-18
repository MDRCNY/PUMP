# Run pump_power on varying values of parameters (grid function)

This extension of \`pump_power()\` will take lists of parameter values
and run \`pump_power()\` on all combinations of these values.

It can only assume the same MDES value for all outcomes due to this.
(I.e., a vector of MDES values will be interpreted as a sequence of
calls to pump_power, one for each MDES value given).

Each parameter in the parameter list can be a list, not scalar. It will
cross all combinations of the list.

## Usage

``` r
pump_power_grid(
  d_m,
  MTP = NULL,
  MDES,
  M = 1,
  nbar,
  J = 1,
  K = 1,
  propZero = NULL,
  numZero = NULL,
  Tbar,
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
  long.table = FALSE,
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

- MDES:

  vector of numeric; This is \*not\* a list of MDES for each outcome,
  but rather a list of MDES to explore. Each value will be assumed held
  constant across all M outcomes.

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

- propZero:

  Proportion of outcomes that have 0 impact (this will be used to
  override numZero, only one can be defined)

- numZero:

  scalar; additional number of outcomes assumed to be zero. Please
  provide NumZero + length(MDES) = M, if length(MDES) is not 1.

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

- long.table:

  TRUE for table with power as rows, correction as columns, and with
  more verbose names. See \`transpose_power_table\`.

- verbose:

  logical; TRUE means print out some text as calls processed. FALSE do
  not.

- drop.unique.columns:

  logical; drop all parameter columns that did not vary across the grid.

- ...:

  extra arguments passed to the underlying pump_power, pump_sample, or
  pump_mdes functions.

## Value

a pumpgridresult object containing power results.

## See also

Other grid functions:
[`pump_mdes_grid()`](https://mdrcny.github.io/PUMP/reference/pump_mdes_grid.md),
[`pump_sample_grid()`](https://mdrcny.github.io/PUMP/reference/pump_sample_grid.md)

## Examples

``` r
g <- pump_power_grid( d_m = "d3.2_m3ff2rc", MTP = c( "HO", "BF" ),
 MDES = 0.10, J = seq(5, 10, 1), M = 5, K = 7, nbar = 58,
 Tbar = 0.50, alpha = 0.15, numCovar.1 = 1,
 numCovar.2 = 1, R2.1 = 0.1, R2.2 = 0.7,
 ICC.2 = 0.25, ICC.3 = 0.25, rho = 0.4, tnum = 1000)
```
