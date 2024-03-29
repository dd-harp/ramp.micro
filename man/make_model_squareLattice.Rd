% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/setup_model.R
\name{make_model_squareLattice}
\alias{make_model_squareLattice}
\title{Setup a square lattice model with one point set offset and nested within the other}
\usage{
make_model_squareLattice(N, kFb, kFq, q_outside = TRUE)
}
\arguments{
\item{N}{the size of the big lattice}

\item{kFb}{a blood feeding dispersal kernel}

\item{kFq}{an egg laying dispersal kernel}

\item{q_outside}{a [logical] switch: if true, habitats define the big lattice}
}
\value{
a model as a compound [list]
}
\description{
Setup a square lattice model with one point set offset and nested within the other
}
