# Update a single pump call to a grid call (grid function)

Take a pumpresult and provide lists of parameters to explore various
versions of the initial scenario.

## Usage

``` r
update_grid(x, ...)
```

## Arguments

- x:

  pump result object.

- ...:

  list of parameters to expand into a grid.

## Value

a pumpgridresult object; result of calling corresponding grid.

## Examples

``` r
pp <- pump_power(d_m = "d2.1_m2fc", MTP = "HO",
  nbar = 200, J = 20, MDES = 0.2, M = 3,
  Tbar = 0.50, alpha = 0.05, numCovar.1 = 5,
  R2.1 = 0.1, ICC.2 = 0.05, rho = 0, tnum = 500)

gd <- update_grid( pp, J = c( 10, 20, 30 ) )
```
