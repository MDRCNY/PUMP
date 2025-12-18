# Provides details about supported package features (core function)

List user options: designs and models (d_m), including what parameters
are relevant for each context; multiple testing procedures; types of
power; design and model parameters.

## Usage

``` r
pump_info(
  topic = c("all", "context", "adjustment", "power", "parameters"),
  comment = TRUE
)
```

## Arguments

- topic:

  string; what kind of info. One of: all, context, adjustment, power,
  parameters.

- comment:

  logical; prints out long description of each design and method.

## Value

list; a list of data frames with information about each topic.

## See also

For more detailed information about user choices, see the manuscript
\<doi:10.18637/jss.v108.i06\>, which includes a detailed Technical
Appendix including information about the designs and models and
parameters.
