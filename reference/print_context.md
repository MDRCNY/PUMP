# Print context (design, model, parameter values) of pumpresult or pumpgridresult (result function)

Print out the context (design and model, with parameter values) of given
pump result or pump grid result object. The "\*\*\*" denotes varying
values in the printout.

## Usage

``` r
print_context(
  x,
  insert_results = FALSE,
  insert_control = FALSE,
  include_SE = TRUE,
  ...
)
```

## Arguments

- x:

  A pumpresult object or pumpgridresult object.

- insert_results:

  Include actual results in the printout.

- insert_control:

  Include the optimizer control parameter information.

- include_SE:

  Include standard errors in the printout.

- ...:

  Extra arguments to pass to print.pumpresult.

## Value

No return value; prints results.
