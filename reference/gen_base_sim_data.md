# Generate base simulated multi-level data (simulation function)

Generates simulated data for multi-level RCTs for pump-supported designs
and models for both unobserved potential outcomes. This function does
not generate treatment assignments or observed outcomes–see
gen_sim_data() for that.

This method takes in a list of necessary data-generating parameters,
following the rest of the package.

This function is beyond the main scope of calculating power, and is
instead used for simulating data. For more info on use, see the
simulation vignette.

## Usage

``` r
gen_base_sim_data(
  param.list,
  pump.object = NULL,
  return.as.dataframe = TRUE,
  no.list = TRUE,
  dgp.params = FALSE
)
```

## Arguments

- param.list:

  list; model parameters such as ICC, R2, etc. See simulation vignette
  for details.

- pump.object:

  A pumpresult object.

- return.as.dataframe:

  TRUE means return list of dataframes, one for each outcome. FALSE
  means return components of the covariates, etc., in a list.

- no.list:

  Only relevant if return.as.dataframe=TRUE. no.list=TRUE means if M=1
  return the dataframe, not a list of length 1. FALSE means return a
  list of length 1, even if there is only 1 outcome.

- dgp.params:

  TRUE means param.list is already converted to DGP parameters, FALSE
  means it needs to be converted via \`convert_params()\`.

## Value

list; potential outcomes given control y0, treatment y1, covariates V.k,
X.jk, C.ijk, or list of dataframes if return.as.dataframe = TRUE.

## See also

gen_sim_data
