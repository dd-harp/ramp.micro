# Bionomic parameters set for the BQ model

Bionomic parameters set for the BQ model

## Usage

``` r
setup_bionomics_BQS(
  model,
  opts = list(),
  pB = 0.96,
  pQ = 0.96,
  pS = 0.96,
  psiB = 0.9,
  psiQ = 0.9,
  psiS = 0.9,
  sigb = 0.1,
  sigq = 0.1,
  sigf = 0.1,
  sigL = 0.1,
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

- pS:

  the probability of surviving a sugar feeding bout

- psiB:

  blood feeding success

- psiQ:

  egg laying success

- psiS:

  sugar feeding success

- sigb:

  transition from blood to sugar feeding

- sigq:

  transition from egg laying to sugar feeding

- sigf:

  transition from egg laying to sugar feeding

- sigL:

  transition from emergence to sugar feeding

- ova:

  female eggs laid during a successful bout, per female

## Value

the parameters, a compound \[list\]
