# Point Sets

``` r
library(ramp.micro)
```

## Introduction

Models in `ramp.micro` are developed around point sets describing the
locations of resources. Point sets are passed to functions that set up
models. For these purposes, the point sets are stored in two columns as
a pair of named vectors.

For example:

``` r
x= c(1,0,1,0)
y= c(0,1,1,0)
xy0 <- cbind(x,y) 
xy0
```

    ##      x y
    ## [1,] 1 0
    ## [2,] 0 1
    ## [3,] 1 1
    ## [4,] 0 0

Note that the points are returned as a pair of named vectors – the first
column is named `x` and the next one `y`

To develop code and facilitate methods development, we developed several
utilities for creating point sets.

## Uniform

The function `unif_xy(n, mn, mx)` generates random uniform point sets
within a square bounding box:

- `mn` $< x <$`mx`

- `mn` $< y <$`mx`

This creates a set of 150 points that are uniformly, randomly
distributed in a 10 unit box:

``` r
xy1 <- unif_xy(50, 0, 10)
head(xy1, 4)
```

    ##              x        y
    ## [1,] 0.8075014 4.611865
    ## [2,] 8.3433304 3.152418
    ## [3,] 6.0076089 1.746759
    ## [4,] 1.5720844 5.315735

Note that `xy1` is $150 \times 2$

``` r
dim(xy1)
```

    ## [1] 50  2

``` r
par(mar = c(2,2,1,1), bty = "o")
plot(xy1)
```

![](point_sets_files/figure-html/unnamed-chunk-5-1.png)

## Lattice

The function `lattice(n, mn, mx)` generates an $n \times n$ lattice
inside a box:

- `mn` $< x <$`mx`

- `mn` $< y <$`mx`

``` r
p2 <- lattice(7, 0, 10)
par(mar = c(2,2,1,1), bty = "o")
plot(p2)
```

![](point_sets_files/figure-html/unnamed-chunk-6-1.png)

``` r
dim(p2)
```

    ## [1] 49  2

## Clusters

Two functions were developed to generate clusters. The first one,
`clusters_xy(xy, nc, vr)` generates a poisson number of points (mean =
`nc`) around a set of seeds. The distance from the center of each seed
is a normal variate with variance `vr.`

``` r
seeds <- unif_xy(23, 0, 50)
p3 <- clusters_xy(seeds, 8, 1)
```

``` r
par(mar = c(2,2,1,1), bty = "o")
plot(p3)
```

![](point_sets_files/figure-html/unnamed-chunk-9-1.png)
