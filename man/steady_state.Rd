% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate.R
\name{steady_state}
\alias{steady_state}
\title{Run the model until it has reached a steady state}
\usage{
steady_state(model, burn = 500, Tx = 50, tol = 0.001)
}
\arguments{
\item{model}{a compound [list] defining a model}

\item{burn}{initial burn time}

\item{Tx}{a run time chunk}

\item{tol}{tolerance for total sum of squared differences after a runtime chunk}
}
\value{
model
}
\description{
Run the model until it has reached a steady state
}
