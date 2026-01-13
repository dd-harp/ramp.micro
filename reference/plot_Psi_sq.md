# Visualize the one-bout dispersal matrices for a BQ model

Visualize the one-bout dispersal matrices for a BQ model

## Usage

``` r
plot_Psi_sq(
  b,
  q,
  s,
  Psi_sq,
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

- Psi_sq:

  one bout dispersal matrix to s from q

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
