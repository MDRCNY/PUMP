# Examine a power curve (result function)

This will give a plot of power vs. MDES or sample size. It can be useful
to see how quickly power changes as a function of these design
parameters. Can be useful to diagnose relatively flat power curves,
where power changes little as a function of MDES or sample size, and can
also be useful to gauge where convergence went poorly.

## Usage

``` r
plot_power_curve(
  pwr,
  plot.points = TRUE,
  all = TRUE,
  low = NULL,
  high = NULL,
  grid.size = 5,
  tnum = 2000,
  breaks = grid.size,
  fit = NULL
)
```

## Arguments

- pwr:

  pumpresult object or data.frame; result from calling pump_sample or
  pump_mdes (or data frame from, e.g., power_curve()).

- plot.points:

  logical; whether to plot individually tested points on curve.

- all:

  logical; if TRUE, merge in the search path from the original search.

- low:

  scalar; low range for curve.

- high:

  scalar; high range for the curve.

- grid.size:

  scalar; number of points to calculate power for.

- tnum:

  scalar; number of iterations to calculate power at each grid point.

- breaks:

  scalar; the desired number of tick marks on the axes.

- fit:

  a four parameter bounded logistic curve (if NULL will fit one to
  passed points).

## Value

plot; a ggplot object of power across values.
