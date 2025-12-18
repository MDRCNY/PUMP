# Validates user inputs

This functions takes in a list of user inputs. Depending on the inputs,
it produces errors or warnings, and at times modifies inputs if
necessary.

## Usage

``` r
validate_inputs(
  d_m,
  params.list,
  power.call = FALSE,
  mdes.call = FALSE,
  ss.call = FALSE,
  verbose = TRUE,
  multi.MTP.ok = FALSE
)
```

## Arguments

- d_m:

  a single RCT d_m (see list/naming convention)

- params.list:

  a list of parameters input by a user

- power.call:

  flag for power estimation

- mdes.call:

  flag for MDES estimation

- ss.call:

  flag for sample size estimation

- verbose:

  whether to print out warnings

- multi.MTP.ok:

  whether validation allows for multiple MTPs

## Value

params.list
