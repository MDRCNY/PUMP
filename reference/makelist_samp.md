# Convert multi-outcome data structure to dataframe for each outcome.

Given the simulated multi-outcome structure, make a list of complete
(tidy) rectangular datasets, one for each outcome.

## Usage

``` r
makelist_samp(samp.obs, T.x = NULL, include_POs = FALSE)
```

## Arguments

- samp.obs:

  a single iteration of observed data

- T.x:

  vector of treatment assignments

- include_POs:

  Include columns for the potential outcomes in addition to the observed
  outcome.

## Value

List of dataframes.
