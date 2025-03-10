% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_VC.R
\name{plot_dispersal_VV}
\alias{plot_dispersal_VV}
\title{Visualize lifetime transmission from a site by mosquito populations}
\usage{
plot_dispersal_VV(
  model,
  cx_b = 2.5,
  cx_q = 0.3,
  min_edge_frac = 0.01,
  r = 0.01,
  arw_lng = 0.05,
  lwd = 2,
  lamp = 1,
  arw_clr = "springgreen4",
  seg_clr = "firebrick3"
)
}
\arguments{
\item{model}{a compound [list] defining a model}

\item{cx_b}{the maximum cex for blood feeding sites}

\item{cx_q}{the maximum cex for egg laying sites}

\item{min_edge_frac}{the fraction of the mass to plot}

\item{r}{the radius of a ring around destination points}

\item{arw_lng}{the arrow length}

\item{lwd}{scale the line width}

\item{lamp}{arrow width scaling factor}

\item{arw_clr}{the color to draw the arrow (asymetric part)}

\item{seg_clr}{the color to draw the segment (symmetric part)}
}
\value{
invisible(NULL)
}
\description{
Visualize lifetime transmission from a site by mosquito populations
}
