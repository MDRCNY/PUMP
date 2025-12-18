# Plot a grid pump mdes object

Plot a grid pump mdes object

## Usage

``` r
# S3 method for class 'pumpgridresult.mdes'
plot(
  x,
  power.definition = NULL,
  var.vary,
  color = "MTP",
  lines = TRUE,
  include.title = FALSE,
  ...
)
```

## Arguments

- x:

  pumpgridresult object.

- power.definition:

  string; definition of power to plot. If NULL, plot all definitions as
  a facet wrap.

- var.vary:

  string; variable to vary on X axis. If NULL, and only one thing
  varies, then it will default to single varying parameter.

- color:

  string; Group lines by this element to make an interaction plot
  (default "MTP", giving one curve for each MTP).

- lines:

  logical; TRUE means connect dots with lines on the plots. FALSE means
  no lines.

- include.title:

  logical; whether to include/exclude title (if planning a facet wrap,
  for example).

- ...:

  additional parameters.
