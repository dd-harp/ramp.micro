# Dispersal parameters for the BQS model

Dispersal parameters for the BQS model

## Usage

``` r
setup_dispersal_BQS(
  model,
  opts = list(),
  kFb,
  kFq,
  kFs,
  wb = 1,
  wq = 1,
  ws = 1,
  stayB = 0.1,
  stayQ = 0.1,
  stayS = 0.1
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

- kFs:

  a function that weights points by distance for sugar searching

- wb:

  search weights for blood feeding sites

- wq:

  search weights for habitats

- ws:

  search weights for sugar sites

- stayB:

  the fraction of blood feeding bouts that stay

- stayQ:

  the fraction of egg laying bouts that stay

- stayS:

  the fraction of sugar feeding bouts that stay

## Value

the model, a compound \[list\]
