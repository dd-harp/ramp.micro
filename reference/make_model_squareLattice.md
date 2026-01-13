# Setup a square lattice model with one point set offset and nested within the other

Setup a square lattice model with one point set offset and nested within
the other

## Usage

``` r
make_model_squareLattice(N, kFb, kFq, q_outside = TRUE)
```

## Arguments

- N:

  the size of the big lattice

- kFb:

  a blood feeding dispersal kernel

- kFq:

  an egg laying dispersal kernel

- q_outside:

  a \[logical\] switch: if true, habitats define the big lattice

## Value

a model as a compound \[list\]
