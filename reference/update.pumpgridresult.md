# Update a pump grid call, tweaking some parameters (core function)

Works on objects returned by \`update_grid()\`; calls \`update_grid()\`.

## Usage

``` r
# S3 method for class 'pumpgridresult'
update(object, ...)
```

## Arguments

- object:

  A pumpgridresult object.

- ...:

  Additional arguments, i.e., the arguments you would pass to the
  \`pump_power()\`, \`pump_mdes()\` and \`pump_sample()\`, that will
  replace the existing parameters of the object.

## See also

\[update_grid()\]
