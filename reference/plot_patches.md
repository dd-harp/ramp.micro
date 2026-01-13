# Visualize the patchespopulation approximation

Visualize the patchespopulation approximation

## Usage

``` r
plot_patches(
  model,
  i,
  cut = NULL,
  f_color = viridis::turbo,
  stretch = 0.1,
  lwd = 2,
  bbend = 3,
  mtl = NULL
)
```

## Arguments

- model:

  a ramp.micro model object

- i:

  the graph

- cut:

  a cutat argument, the \# of communities

- f_color:

  a function that returns a list of colors (e.g. viridis::turbo)

- stretch:

  make the hull slightly larger or slightly smaller

- lwd:

  a plotting argument

- bbend:

  the bend in the arrows

- mtl:

  a title for the graph

## Value

a ramp.micro model object
