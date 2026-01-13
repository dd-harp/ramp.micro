# Plot lifetime egg dispersal, per mosquito

Plot lifetime egg dispersal, per mosquito

## Usage

``` r
plot_dispersal_G(
  model,
  cx_b = 0.3,
  cx_q = 2,
  min_edge_frac = 0.01,
  r = 0.01,
  arw_lng = 0.05,
  lwd = 2,
  lamp = 1,
  seg_clr = "lightblue",
  arw_clr = "salmon"
)
```

## Arguments

- model:

  a model defined as a compound \[list\]

- cx_b:

  the maximum cex for blood feeding sites

- cx_q:

  the maximum cex for egg laying sites

- min_edge_frac:

  the fraction of the mass to plot

- r:

  the radius of a ring around destination points

- arw_lng:

  the arrow length

- lwd:

  scale the line width

- lamp:

  arrow width scaling factor

- seg_clr:

  the color to draw the segment (symmetric part)

- arw_clr:

  the color to draw the arrow (asymetric part)

## Value

invisible(NULL)
