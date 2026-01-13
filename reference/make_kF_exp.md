# Make an exponential function for weight by distance

This returns a function of the form \$\$F_w (d, \omega=1) = \omega_j
e^{-k \left( \frac{d\_{i,j}}{s}\right)^\gamma}\$\$ where \\s\\ and
\\\gamma\\ are shape parameters, \\k\\ is the rate parameter, and
\\\omega\\ is a weight.

In effect, \\s\\ is the location of a shoulder, and for \\\gamma\>1\\,
the decay is slower for \\d\<s\\.

The function returned accepts \\\omega\\ as an optional argument so that
it can be passed at the time of simulation.

By default, the function returns scaled values – the maximum is 1.

## Usage

``` r
make_kF_exp(k = 1, s = 2, gamma = 1)
```

## Arguments

- k:

  decay by distance

- s:

  a scale parameter

- gamma:

  a shape parameter

## Value

a function

## Examples

``` r
kF1 = make_kF_exp(k=1, s=1, gamma=1.5)
kF2 = make_kF_exp(k=2, s=0.1, gamma=2)
dd = seq(0, 2, by = 0.01)
plot(dd, kF1(dd), type = "l", ylab = "Weight", xlab = "Distance")
lines(dd, kF2(dd))
```
