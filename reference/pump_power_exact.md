# Calculate power theoretically for M=1 situations

Calculate power theoretically for M=1 situations

## Usage

``` r
pump_power_exact(MDES, SE, df, alpha, two.tailed)
```

## Arguments

- MDES:

  MDES (single number)

- SE:

  Calculated SE of the estimator

- df:

  Degrees of freedom of the estimator (often approximated).

- alpha:

  Alpha for the planned test.

- two.tailed:

  TRUE/FALSE Two- or one-sided test?

## Value

Single row Tibble with columns or power, SE, and DF. MTP column with
value of "None".
