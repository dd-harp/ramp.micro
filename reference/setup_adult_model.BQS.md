# Setup a BQS model for adult mosquitoes

Setup a BQS model for adult mosquitoes

## Usage

``` r
# S3 method for class 'BQS'
setup_adult_model(
  model,
  b,
  q,
  s,
  dispersal_opts = list(),
  bionomic_opts = list(),
  eip = 15
)
```

## Arguments

- model:

  a model defined as a compound \[list\]

- b:

  a point set defining blood feeding sites

- q:

  a point set defining egg laying sites

- s:

  a point set defining sugar feeding sites

- dispersal_opts:

  a \[list\] to overwrite defaults

- bionomic_opts:

  a \[list\] to overwrite defaults

- eip:

  the extrinsic incubation period

## Value

a \[list\] defining a BQ-class adult model
