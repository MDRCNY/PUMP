# Calculate degrees of freedom (support function)

Given sample sizes, return the used degrees of freedom (frequently
conservative) for the design and model.

## Usage

``` r
calc_df(d_m, J, K, nbar, numCovar.1, numCovar.2, numCovar.3, validate = TRUE)
```

## Arguments

- d_m:

  string; a single context, which is a design and model code. See
  pump_info() for list of choices.

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

- numCovar.1:

  scalar; number of level 1 (individual) covariates.

- numCovar.2:

  scalar; number of level 2 (school) covariates.

- numCovar.3:

  scalar; number of level 3 (district) covariates.

- validate:

  logical; whether or not to validate if output df is \<= 0.

## Value

scalar; degrees of freedom for the context.
