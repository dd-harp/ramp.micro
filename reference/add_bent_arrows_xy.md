# Draw bent arrows from one point set to another

Draw bent arrows from one point set to another

## Usage

``` r
add_bent_arrows_xy(
  xy_launch,
  xy_land,
  M,
  mnwd = 0.05,
  bbend = 0,
  endd = 0.6,
  adj = 2,
  clr = "red"
)
```

## Arguments

- xy_launch:

  origin point set

- xy_land:

  destination point set

- M:

  the dispersal matrix

- mnwd:

  the minimum width to plot

- bbend:

  a parameter to bend arrows

- endd:

  relative position of arrowhead

- adj:

  set the maximum cex for points

- clr:

  the color

## Value

invisible(NULL)
