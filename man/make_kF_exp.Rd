% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/setup_kernels.R
\name{make_kF_exp}
\alias{make_kF_exp}
\title{Make an exponential function for weight by distance}
\usage{
make_kF_exp(k = 1, s = 2, gamma = 1)
}
\arguments{
\item{k}{decay by distance}

\item{s}{a scale parameter}

\item{gamma}{a shape parameter}
}
\value{
a function
}
\description{
This returns a function of the form
\deqn{F_w (d, \omega=1) = \omega_j e^{-k \left( \frac{d_{i,j}}{s}\right)^\gamma}}
where \eqn{s} and \eqn{\gamma} are shape parameters, \eqn{k} is the rate
parameter, and \eqn{\omega} is a weight.

In effect, \eqn{s} is the location of a shoulder, and for \eqn{\gamma>1}, the decay is
slower for \eqn{d<s}.

The function returned accepts \eqn{\omega} as
an optional argument so that it can be passed at the time of simulation.

By default, the function returns scaled values -- the maximum is 1.
}
\examples{
kF1 = make_kF_exp(k=1, s=1, gamma=1.5)
kF2 = make_kF_exp(k=2, s=0.1, gamma=2)
dd = seq(0, 2, by = 0.01)
plot(dd, kF1(dd), type = "l", ylab = "Weight", xlab = "Distance")
lines(dd, kF2(dd))
}
