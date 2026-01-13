# Visualize a dispersal matrix among two point sets

Visualize a dispersal matrix among two point sets

## Usage

``` r
plot_matrix_xy(
  xy_launch,
  xy_land,
  M,
  mx_pt_sz1 = 0.7,
  mx_pt_sz2 = 1.5,
  pt_clr1 = "salmon",
  pt_clr2 = "lightblue",
  min_edge_frac = 0.01,
  r = 0,
  arw_lng = 0.02,
  arw_lwd = 2,
  lamp = 1,
  arw_clr = "tomato",
  mtl = ""
)
```

## Arguments

- xy_launch:

  origin point set

- xy_land:

  destination point set

- M:

  the dispersal matrix

- mx_pt_sz1:

  the maximum cex for the launch point sest

- mx_pt_sz2:

  the maximum cex for the landing point set

- pt_clr1:

  point color for launch

- pt_clr2:

  point color for landing

- min_edge_frac:

  the minimum fraction to plot edge

- r:

  the radius of a ring around destination points

- arw_lng:

  the arrow length

- arw_lwd:

  scale the line width

- lamp:

  arrow width scaling factor

- arw_clr:

  the arrow color

- mtl:

  = plot title

## Value

no visible return value
