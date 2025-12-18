# Convert power table from wide to long (result function)

Transform table returned from pump_power to a long format table or to a
wide format table.

## Usage

``` r
transpose_power_table(power_table, M = NULL)
```

## Arguments

- power_table:

  pumpresult object for a power result (not mdes or sample). (It can
  also take a raw dataframe of the wide table to convert to long, as an
  internal helper method.)

- M:

  scalar; set if power_table is a data.frame without set number of
  outcomes. Usually ignore this.

## Value

data.frame of power results in long format.
