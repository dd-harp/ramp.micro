# Dispersal parameters for the BQ model

Dispersal parameters for the BQ model

## Usage

``` r
setup_dispersal_BQ(
  model,
  opts = list(),
  kFb,
  kFq,
  wb = 1,
  wq = 1,
  stayB = 0.1,
  stayQ = 0.1
)
```

## Arguments

- model:

  a model defined as a compound \[list\]

- opts:

  a \[list\] of values that overwrite the defaults

- kFb:

  a function that weights points by distance for blood searching

- kFq:

  a function that weights points by distance for habitat searching

- wb:

  search weights for blood feeding sites

- wq:

  search weights for habitats

- stayB:

  the fraction of blood feeding bouts that stay

- stayQ:

  the fraction of egg laying bouts that stay

## Value

the model, a compound \[list\]
