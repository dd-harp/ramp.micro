# Visualize the one-bout dispersal matrices for a BQ model

Visualize the one-bout dispersal matrices for a BQ model

## Usage

``` r
plot_all_Psi_BQ(
  b,
  q,
  Psi_bb,
  Psi_qb,
  Psi_bq,
  Psi_qq,
  cx_D = 2,
  cx_S = 0.3,
  min_edge_frac = 0.01,
  r = 0.01,
  arw_lng = 0.05,
  lwd = 2
)
```

## Arguments

- b:

  blood feeding sites point set

- q:

  egg laying sites point set

- Psi_bb:

  one bout dispersal matrix from b to b

- Psi_qb:

  one bout dispersal matrix from q to b

- Psi_bq:

  one bout dispersal matrix from b to q

- Psi_qq:

  one bout dispersal matrix from q to q

- cx_D:

  the maximum cex for the destination

- cx_S:

  the maximum cex for the source

- min_edge_frac:

  the fraction of the mass to plot

- r:

  the radius of a ring around destination points

- arw_lng:

  the arrow length

- lwd:

  scale the line width

## Value

no visible return value
