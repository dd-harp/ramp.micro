# Visualize a dispersal matrix from one point set to itself

Visualize a dispersal matrix from one point set to itself

## Usage

``` r
plot_matrix_xx(
  xy,
  M,
  mx_pt_sz = 2,
  pt_clr = "tomato",
  min_edge_frac = 0.01,
  r = 0.01,
  arw_lng = 0.05,
  lwd = 1,
  lamp = 1,
  arw_clr = "#e2739655",
  seg_clr = "#00000022",
  mtl = ""
)
```

## Arguments

- xy:

  the point set

- M:

  the dispersal matrix

- mx_pt_sz:

  set cex for the largest point

- pt_clr:

  the point color

- min_edge_frac:

  the minimum fraction to plot edge

- r:

  the radius of a ring around destination points

- arw_lng:

  the arrow length

- lwd:

  set the maximum line width

- lamp:

  arrow width scaling factor

- arw_clr:

  the color to draw the arrow (asymetric part)

- seg_clr:

  the color to draw the segment (symmetric part)

- mtl:

  = plot title

## Value

no visible return value
