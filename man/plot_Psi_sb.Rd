% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/viz_Psi_matrices.R
\name{plot_Psi_sb}
\alias{plot_Psi_sb}
\title{Visualize the one-bout dispersal matrices for a BQ model}
\usage{
plot_Psi_sb(
  b,
  q,
  s,
  Psi_sb,
  cx_D = 2,
  cx_S = 0.3,
  min_edge_frac = 0.01,
  r = 0.01,
  arw_lng = 0.05,
  lwd = 2
)
}
\arguments{
\item{b}{blood feeding sites point set}

\item{q}{egg laying sites point set}

\item{s}{sugar feeding sites point set}

\item{Psi_sb}{one bout dispersal matrix to s from b}

\item{cx_D}{the maximum cex for the destination}

\item{cx_S}{the maximum cex for the source}

\item{min_edge_frac}{the fraction of the mass to plot}

\item{r}{the radius of a ring around destination points}

\item{arw_lng}{the arrow length}

\item{lwd}{scale the line width}
}
\value{
no visible return value
}
\description{
Visualize the one-bout dispersal matrices for a BQ model
}
