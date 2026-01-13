# Setup a behavioral state, micro-simulation model

Setup a behavioral state, micro-simulation model

## Usage

``` r
setup_model(
  b,
  q,
  s = c(),
  kFb = NULL,
  kFq = NULL,
  kFs = NULL,
  Mname = "BQ",
  Lname = "basicL",
  dispersal_opts = list(),
  bionomic_opts = list(),
  aquatic_opts = list(),
  M0_opts = list(),
  L0_opts = list(),
  eip = 15
)
```

## Arguments

- b:

  a point set defining blood feeding sites

- q:

  a point set defining egg laying sites

- s:

  a point set defining sugar feeding sites, with a default null value

- kFb:

  a kernel shape for blood searching

- kFq:

  a kernal shape for aquatic habitat searching

- kFs:

  a kernel shape for sugar site searching

- Mname:

  the adult model name

- Lname:

  the aquatic model name

- dispersal_opts:

  a \[list\] to

- bionomic_opts:

  a \[list\] to overwrite defaults

- aquatic_opts:

  a \[list\] to overwrite defaults

- M0_opts:

  options to overwrite defaults

- L0_opts:

  options to overwrite defaults

- eip:

  the extrinsic incubation period

## Value

a \[list\] defining a BQ-class adult model
