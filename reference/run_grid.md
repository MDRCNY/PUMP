# Run grid across any of the core pump functions

Run grid across any of the core pump functions

## Usage

``` r
run_grid(
  args,
  pum_function,
  verbose = FALSE,
  drop.unique.columns,
  ...,
  use.furrr = FALSE
)
```

## Arguments

- args:

  list of scenario arguments

- pum_function:

  pump_mdes, pump_sample, pump_power

- verbose:

  print out detailed diagnostics

- drop.unique.columns:

  logical; drop all parameter columns that did not vary across the grid.

- ...:

  Extra arguments passed to the underlying pump_power, pump_sample, or
  pump_mdes functions.

- use.furrr:

  not currently implemented; whether to use furr package for
  parallelization
