# Visualize lifetime transmission from a site by mosquito populations

Visualize lifetime transmission from a site by mosquito populations

## Usage

``` r
plot_dispersal_VV(
  model,
  cx_b = 2.5,
  cx_q = 0.3,
  min_edge_frac = 0.01,
  r = 0.01,
  arw_lng = 0.05,
  lwd = 2,
  lamp = 1,
  arw_clr = "springgreen4",
  seg_clr = "firebrick3"
)
```

## Arguments

- model:

  a compound \[list\] defining a model

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

- arw_clr:

  the color to draw the arrow (asymetric part)

- seg_clr:

  the color to draw the segment (symmetric part)

## Value

invisible(NULL)
