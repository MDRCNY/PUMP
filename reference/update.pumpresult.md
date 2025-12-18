# Update a pump call, tweaking some parameters (core function)

Works on objects returned by pump_power(), pump_mdes(), or
pump_sample(). One of the optional parameters can be a \`type =
something\` argument, where the "something" is either "power", "sample",
or "mdes", if the call should be shifted to a different pump call
(pump_power, pump_sample, or pump_mdes, respectively).

## Usage

``` r
# S3 method for class 'pumpresult'
update(object, type = NULL, ...)
```

## Arguments

- object:

  pump result object.

- type:

  string; can be "power", "mdes" or "sample", sets the type of the
  updated call (can be different from original).

- ...:

  parameters as specified in \`pump_power\`, \`pump_mdes\`, and
  \`pump_sample\` that should be overwritten.

## Value

a pumpresult object: results of a new call using parameters of old
object with newly specified parameters replaced.

## Examples

``` r
ss <- pump_sample( d_m = "d2.1_m2fc", MTP = "HO",
  typesample = "J", nbar = 200, power.definition = "min1",
  M = 5, MDES = 0.05, target.power = 0.5, tol = 0.05,
  Tbar = 0.50, alpha = 0.05, numCovar.1 = 5, R2.1 = 0.1,
  ICC.2 = 0.15, rho = 0, final.tnum = 1000 )

up <- update(ss, nbar = 40, tnum = 2000 )
```
