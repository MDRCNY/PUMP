# Determine next point to check for correct power level.

Extract roots from power curve fit to given evaluated points

## Usage

``` r
find_best(test.pts, target.power, gamma = 1.5)
```

## Arguments

- test.pts:

  power evaluated at different points

- target.power:

  goal power

- gamma:

  Number \> 1. The amount we can extend our search range up (gamma) or
  down (1/gamma) if we want to search outside our current boundaries. 1
  means do not extend (which will potentially cause search to stall, not
  recommended).

## Value

List of estimate of when curve reaches target.power, derivative of curve
at that point, and parameters of the fit curve.
