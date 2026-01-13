# Set up the aquatic population model of class basicL

Set up the aquatic population model of class basicL

## Usage

``` r
setup_aquatic_model_basicL(
  model,
  opts = list(),
  pL = 0.9,
  zeta = 0.01,
  theta = 0.9,
  xi = 0
)
```

## Arguments

- model:

  a model defined as a compound \[list\]

- opts:

  a \[list\] of values that overwrite the defaults

- pL:

  base survival fraction

- zeta:

  reduced survival in response to mean crowding

- theta:

  base maturation fraction

- xi:

  delayed maturation in response to mean crowding

## Value

the model, a compound \[list\]
