# Calculate a power curve for sample size or mdes.

For a grid of points based on a passed sample or mdes pumpresult,
estimate power.

## Usage

``` r
estimate_power_curve(
  p,
  low = NULL,
  high = NULL,
  high.max = 1000,
  grid.size = 5,
  tnum = 2000
)
```

## Arguments

- p:

  pumpresult object

- low:

  Low end of grid

- high:

  High end of grid

- high.max:

  If no high provided, maximum possible high

- grid.size:

  Number of points in grid

- tnum:

  the number of test statistics (samples)

## Value

List of powers for grid of points.
