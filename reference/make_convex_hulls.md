# Make a set of convex hulls to visualize the community

Make a set of convex hulls to visualize the community

## Usage

``` r
make_convex_hulls(
  model,
  graph,
  cut = NULL,
  clrs = NULL,
  f_color = viridis::turbo
)
```

## Arguments

- model:

  a model defined as a compound \[list\]

- graph:

  a graphwork object

- cut:

  if !null, the number of communities for igraph::cut_at

- clrs:

  a set of colors

- f_color:

  a function that returns a list of colors (e.g. viridis::turbo)

## Value

the graphwork object with convex hulls attached
