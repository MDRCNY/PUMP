# Generate treatment assignment vector (simulation function)

Given a RCT design and supporting information, generates treatment
assignments for each student.

This function is beyond the main scope of calculating power, and is
instead used for simulating data. For more info on use, see the
simulation vignette.

## Usage

``` r
gen_T.x(d_m, S.id, D.id, Tbar)
```

## Arguments

- d_m:

  string; design and model.

- S.id:

  vector; school assignments.

- D.id:

  vector; district assignments.

- Tbar:

  scalar; probability of treatment assignment.

## Value

vector; treatment assignments for each unit.
