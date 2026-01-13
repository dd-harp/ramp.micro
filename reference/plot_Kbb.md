# Plot the matrix \\K\_{b \leftarrow b}\\

Plot the matrix \\K\_{b \leftarrow b}\\

## Usage

``` r
plot_Kbb(
  model,
  cx_b = 2,
  cx_q = 0.3,
  min_edge_frac = 0.01,
  r = 0.02,
  arw_lng = 0.002,
  lwd = 2,
  arw_clr = "#e2739655",
  seg_clr = "#00000022",
  clr_b = "#cc444bCC",
  clr_q = "#4361eeCC"
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

  the fraction of the mass to plot

- r:

  the radius of a ring around destination points

- arw_lng:

  the arrow length

- lwd:

  scale the line width

- arw_clr:

  the color to draw the arrow (asymetric part)

- seg_clr:

  the color to draw the segment (symmetric part)

- clr_b:

  color for blood feeding sites

- clr_q:

  color for egg laying sites
