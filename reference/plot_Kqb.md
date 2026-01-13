# Plot the matrix \\K\_{q \leftarrow b}\\ dispersal from \\\left\\b\right\\\\ to \\\left\\q\right\\\\

Plot the matrix \\K\_{q \leftarrow b}\\ dispersal from
\\\left\\b\right\\\\ to \\\left\\q\right\\\\

## Usage

``` r
plot_Kqb(
  model,
  cx_b = 0.3,
  cx_q = 2,
  min_edge_frac = 0.01,
  r = 0.02,
  arw_lng = 0.002,
  lwd = 2,
  clr_K = "#4361eeCC",
  clr_b = "red",
  clr_q = "darkblue"
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

- clr_K:

  the color for Kqb arrows

- clr_b:

  the color for b points

- clr_q:

  the color for q points

## Value

invisible(NULL)
