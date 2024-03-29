% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adult-BQS.R
\name{setup_dispersal_BQS}
\alias{setup_dispersal_BQS}
\title{Dispersal parameters for the BQS model}
\usage{
setup_dispersal_BQS(
  model,
  opts = list(),
  kFb,
  kFq,
  kFs,
  wb = 1,
  wq = 1,
  ws = 1,
  stayB = 0.1,
  stayQ = 0.1,
  stayS = 0.1
)
}
\arguments{
\item{model}{a model defined as a compound [list]}

\item{opts}{a [list] of values that overwrite the defaults}

\item{kFb}{a function that weights points by distance for blood searching}

\item{kFq}{a function that weights points by distance for habitat searching}

\item{kFs}{a function that weights points by distance for sugar searching}

\item{wb}{search weights for blood feeding sites}

\item{wq}{search weights for habitats}

\item{ws}{search weights for sugar sites}

\item{stayB}{the fraction of blood feeding bouts that stay}

\item{stayQ}{the fraction of egg laying bouts that stay}

\item{stayS}{the fraction of sugar feeding bouts that stay}
}
\value{
the model, a compound [list]
}
\description{
Dispersal parameters for the BQS model
}
