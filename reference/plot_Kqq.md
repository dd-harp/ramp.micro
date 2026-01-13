# Plot the matrix \\K\_{q \leftarrow q}\\

Plot the matrix \\K\_{q \leftarrow q}\\

## Usage

``` r
plot_Kqq(
  model,
  cx_b = 0.03,
  cx_q = 2,
  min_edge_frac = 0.01,
  r = 0.02,
  arw_lng = 0.002,
  lwd = 2,
  arw_clr = "#4361eeCC",
  clr_qA = "#abc4ff55",
  seg_clr = "#00000022",
  clr_qB = "#858ae399",
  clr_b = "#cc444bCC"
)
```

## Arguments

- model:

  a model defined as a compound \[list\]

- cx_b:

  set the maximum cex for b points

- cx_q:

  set the maximum cex for q points

- min_edge_frac:

  the minimum fraction to plot edge

- r:

  the radius of a ring around destination points

- arw_lng:

  the arrow length

- lwd:

  scale the line width

- arw_clr:

  the color to draw the arrow (asymetric part)

- clr_qA:

  color for egg laying sites

- seg_clr:

  the color to draw the segment (symmetric part)

- clr_qB:

  color for egg laying sites

- clr_b:

  color for blood feeding sites

## Value

invisible(NULL)
