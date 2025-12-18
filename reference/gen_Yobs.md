# Generate observed outcomes (simulation function)

Takes in a full dataset of both observed and latent potential outcomes
and the treatment assignment vector, and returns only the observed
outcomes.

This function is beyond the main scope of calculating power, and is
instead used for simulating data. For more info on use, see the
simulation vignette.

## Usage

``` r
gen_Yobs(full.data, T.x)
```

## Arguments

- full.data:

  data.frame; full dataset of potential outcomes.

- T.x:

  vector; binary assignment to treat/control.

## Value

vector; observed outcomes
