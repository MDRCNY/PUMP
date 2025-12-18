# Examine search path of a power search (result function)

This will give triple-plots about how the search narrowed down into the
final estimate. Can be useful to gauge where convergence went poorly.

## Usage

``` r
plot_power_search(pwr, fit = NULL, target.line = NULL)
```

## Arguments

- pwr:

  pumpresult object; result from a pump_sample or pump_mdes call.

- fit:

  a fitted curve to the search.

- target.line:

  scalar; if non-NULL, add a reference line for the true power (if
  known, e.g., from a pump_power call).

## Value

plot; a ggplot object (a ggpubr arrangement of 3 plots, technically) of
the search path.
