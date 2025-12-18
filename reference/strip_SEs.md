# Remove SE and df columns from (wide) power table

This is used to reduce the info on a power table before pivoting to long
format.

## Usage

``` r
strip_SEs(power_table)
```

## Arguments

- power_table:

  Dataframe (power result object).

## Value

Changed dataframe with all columns starting with SE or df dropped.
