# Plot the ouptuts of a graph using

Plot the ouptuts of a graph using

## Usage

``` r
plot_subgraph(
  model,
  graph,
  cut = NULL,
  alg = "wt",
  f_color = viridis::turbo,
  min_edge_frac = 0.01,
  cx = 2,
  pw = 1,
  mtl = "",
  stretch = 0.1,
  r = 0.02,
  arw_lng = 0.05,
  lwd = 2
)
```

## Arguments

- model:

  a \`ramp.micro\` model object

- cut:

  optional arguent for cut_at

- alg:

  walktrap = "wt" or greedy = "gr"

- f_color:

  a function that returns a list of colors (e.g. viridis::turbo)

- min_edge_frac:

  the fraction of the mass to plot

- cx:

  cex for plotting points

- pw:

  power relationship for scaling point size: pw=1 is linear

- mtl:

  a plot title

- stretch:

  make the hull slightly larger or slightly smaller

- r:

  the radius of a ring around destination points

- arw_lng:

  the arrow length

- lwd:

  scale the line width

- graphs:

  a graphswork object

## Value

invisible(NULL)
