# Bionomic parameters set for the BQ model

Bionomic parameters set for the BQ model

## Usage

``` r
setup_bionomics_BQ(
  model,
  opts = list(),
  pB = 0.96,
  pQ = 0.96,
  psiB = 0.9,
  psiQ = 0.9,
  ova = 20
)
```

## Arguments

- model:

  a model defined as a compound \[list\]

- opts:

  a \[list\] of values that overwrite the defaults

- pB:

  the probability of surviving a blood feeding bout

- pQ:

  the probability of surviving an egg laying bout

- psiB:

  blood feeding success

- psiQ:

  egg laying success

- ova:

  female eggs laid during a successful bout, per female

## Value

the parameters, a compound \[list\]
