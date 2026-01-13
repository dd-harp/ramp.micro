# Visualize the one-bout dispersal matrices for a BQS model

Visualize the one-bout dispersal matrices for a BQS model

## Usage

``` r
plot_all_Psi_BQS(
  b,
  q,
  s,
  Psi_bb,
  Psi_qb,
  Psi_sb,
  Psi_bq,
  Psi_qq,
  Psi_sq,
  Psi_bs,
  Psi_qs,
  Psi_ss,
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

- s:

  sugar feeding sites point set

- Psi_bb:

  one bout dispersal matrix from b to b

- Psi_qb:

  one bout dispersal matrix from q to b

- Psi_sb:

  one bout dispersal matrix from s to b

- Psi_bq:

  one bout dispersal matrix from b to q

- Psi_qq:

  one bout dispersal matrix from q to q

- Psi_sq:

  one bout dispersal matrix from s to q

- Psi_bs:

  one bout dispersal matrix from b to s

- Psi_qs:

  one bout dispersal matrix from q to s

- Psi_ss:

  one bout dispersal matrix from s to s

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
