% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/K_matrices.R
\name{plot_Kbq}
\alias{plot_Kbq}
\title{Plot the matrix Kbq: dispersal from {q} to {b}}
\usage{
plot_Kbq(
  model,
  max_pt_sz = 2,
  min_edge_frac = 0.01,
  r = 0.02,
  arw_lng = 0.002,
  lwd = 2,
  clr_K = "#fe5f55CC",
  clr_b = "darkred",
  clr_q = "#858ae399"
)
}
\arguments{
\item{model}{a model defined as a compound [list]}

\item{max_pt_sz}{set the maximum cex for points}

\item{min_edge_frac}{the fraction of the mass to plot}

\item{r}{the radius of a ring around destination points}

\item{arw_lng}{the arrow length}

\item{lwd}{scale the line width}

\item{clr_K}{the color for Kqb arrows}

\item{clr_b}{the color for b points}

\item{clr_q}{the color for q points}
}
\value{
invisible(NULL)
}
\description{
Plot the matrix Kbq: dispersal from {q} to {b}
}
