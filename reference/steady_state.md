# Run the model until it has reached a steady state

Run the model until it has reached a steady state

## Usage

``` r
steady_state(model, burn = 500, Tx = 50, tol = 0.001)
```

## Arguments

- model:

  a compound \[list\] defining a model

- burn:

  initial burn time

- Tx:

  a run time chunk

- tol:

  tolerance for total sum of squared differences after a runtime chunk

## Value

model
