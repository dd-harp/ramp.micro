# Make graph object

From a matrix, M, return a graph object the clusters from the walktrap
algorithm & for the symmetric graph defined by M + t(M) the clusters for
greedy_clusters

## Usage

``` r
make_graph_obj(M, type = "b", tag = "")
```

## Arguments

- M:

  a matrix

- type:

  "b" for blood feeding; "q" for egg laying

- tag:

  a string for plotting

## Value

a graph object
