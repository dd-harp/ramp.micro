% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_points.R
\name{add_points_ss}
\alias{add_points_ss}
\title{Add sugar feeding sites with and without the self loop}
\usage{
add_points_ss(
  s,
  M,
  pw = 1,
  max_pt_sz = 2,
  colA = "#8bc74266",
  colB = "#8bc742CC"
)
}
\arguments{
\item{s}{a point set defining sugar feeding sites}

\item{M}{a matrix describing movement to s}

\item{pw}{power relationship for scaling point size: pw=1 is linear}

\item{max_pt_sz}{set the maximum cex for points}

\item{colA}{color with the self loop}

\item{colB}{color without the self loop}
}
\value{
invisible(NULL)
}
\description{
Add sugar feeding sites with and without the self loop
}
