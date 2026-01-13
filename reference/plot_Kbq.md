# Plot the matrix \\K\_{b \leftarrow q}\\: dispersal from \\\left\\q\right\\\\ to \\\left\\b\right\\\\

Plot the matrix \\K\_{b \leftarrow q}\\: dispersal from
\\\left\\q\right\\\\ to \\\left\\b\right\\\\

## Usage

``` r
plot_Kbq(
  model,
  cx_b = 2,
  cx_q = 0.3,
  min_edge_frac = 0.01,
  r = 0.02,
  arw_lng = 0.002,
  lwd = 2,
  clr_K = "#fe5f55CC",
  clr_b = "darkred",
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

- clr_K:

  the color for Kqb arrows

- clr_b:

  the color for b points

- clr_q:

  the color for q points

## Value

invisible(NULL)
