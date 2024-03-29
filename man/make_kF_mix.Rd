% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kernels.R
\name{make_kF_mix}
\alias{make_kF_mix}
\title{Make an power function that combines an exponential and power function}
\usage{
make_kF_mix(p = 0.001, k = 1, s = 2, gamma = 1, delta = 1)
}
\arguments{
\item{p}{the weight on the power function}

\item{k}{decay by distance}

\item{s}{a scale parameter}

\item{gamma}{a shape parameter}

\item{delta}{the power on distance}
}
\value{
a function
}
\description{
Make an power function that combines an exponential and power function
}
