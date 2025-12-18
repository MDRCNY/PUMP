# Result object for results of grid power calculations

The pumpgridresult object is an S3 class that holds the results from
\`pump_power_grid()\`, \`pump_sample_grid()\`, and \`pump_mdes_grid()\`.

It has several methods that pull different information from this object,
and some printing methods for getting nicely formatted results.

## Usage

``` r
is.pumpgridresult(x)

# S3 method for class 'pumpgridresult'
print(x, header = TRUE, include_SE = FALSE, ...)

# S3 method for class 'pumpgridresult'
summary(object, include_SE = FALSE, ...)
```

## Arguments

- x:

  a pumpgridresult object (except for is.pumpgridresult, where it is a
  generic object to check).

- header:

  logical; FALSE means skip some header info on the result, just print
  the data.frame of actual results.

- include_SE:

  logical; TRUE means include standard errors and df.

- ...:

  extra options passed to print.pumpgridresult

- object:

  object to summarize.

## Value

is.pumpgridresult: TRUE if object is a pumpgridresult object.

print: No return value; prints results.

summary: No return value; prints results.
