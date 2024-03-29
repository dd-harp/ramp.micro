% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/K_matrices.R
\name{plot_Kbb}
\alias{plot_Kbb}
\title{Plot the matrix Kbb}
\usage{
plot_Kbb(
  model,
  max_pt_sz = 2,
  min_edge_frac = 0.01,
  r = 0.02,
  arw_lng = 0.002,
  lwd = 2,
  arw_clr = "#e2739655",
  seg_clr = "#00000022",
  clr_b = "#cc444bCC",
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

\item{arw_clr}{the color to draw the arrow (asymetric part)}

\item{seg_clr}{the color to draw the segment (symmetric part)}

\item{clr_b}{color for blood feeding sites}

\item{clr_q}{color for egg laying sites}
}
\description{
Plot the matrix Kbb
}
