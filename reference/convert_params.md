# Converts model params into DGP params (simulation function)

Converts user-provided parameters such as ICC and omega into
data-generating parameters for the multilevel random effects model used
to produce simulated data, such as variance values and covariate
coefficients.

This function is beyond the main scope of calculating power, and is
instead used for simulating data. For more info on use, see the
simulation vignette.

## Usage

``` r
convert_params(param.list)
```

## Arguments

- param.list:

  list; model parameters such as ICC, R2, etc.

## Value

list; data-generating parameters.
