
#' Plot lifetime egg dispersal, per mosquito
#'
#' @param model a model defined as a compound [list]
#' @param mx_pt_sz_b the maximum cex for blood feeding sites
#' @param mx_pt_sz_q the maximum cex for egg laying sites
#' @param min_edge_frac the fraction of the mass to plot
#' @param r the radius of a ring around destination points
#' @param arw_lng the arrow length
#' @param lwd scale the line width
#' @param lamp arrow width scaling factor
#' @param arw_clr the color to draw the arrow (asymetric part)
#' @param seg_clr the color to draw the segment (symmetric part)
#'
#' @return invisible(NULL)
#' @export
plot_dispersal_G = function(model,
                    mx_pt_sz_b = 0.7, mx_pt_sz_q = 2,
                    min_edge_frac = 0.01, r=.01, arw_lng=0.05, lwd=2, lamp=1,
                    seg_clr="lightblue", arw_clr="salmon"){
  with(model,{
    par(mar=c(2,2,2,2))
    frame_bq(b, q, mtl = "Lifetime Egg Dispersal, per Mosquito")
    add_points_b(b, max_pt_sz = mx_pt_sz_b)
    add_arrows_xx(q, G, min_edge_frac=min_edge_frac, r=r, arw_lng=arw_lng, lwd=lwd,
                  lamp=lamp, arw_clr=arw_clr, seg_clr=seg_clr)
    add_points_qq(q, G, max_pt_sz = mx_pt_sz_q)
  })
  return(invisible())
}

#' Plot lifetime egg dispersal by a mosquito population
#'
#' @param model a model defined as a compound [list]
#' @param mx_pt_sz_b the maximum cex for blood feeding sites
#' @param mx_pt_sz_q the maximum cex for egg laying sites
#' @param min_edge_frac the fraction of the mass to plot
#' @param r the radius of a ring around destination points
#' @param arw_lng the arrow length
#' @param lwd scale the line width
#' @param lamp arrow width scaling factor
#' @param arw_clr the color to draw the arrow (asymetric part)
#' @param seg_clr the color to draw the segment (symmetric part)
#'
#' @return invisible(NULL)
#' @export
plot_dispersal_GG = function(model,
                    mx_pt_sz_b = 0.7, mx_pt_sz_q = 2,
                    min_edge_frac = 0.01, r=.01, arw_lng=0.05, lwd=2, lamp=1,
                    seg_clr="steelblue", arw_clr="chocolate"){
  with(model,{
    par(mar=c(2,2,2,2))
    frame_bq(b, q, mtl = "Lifetime Egg Dispersal, Population")
    add_points_b(b, max_pt_sz = mx_pt_sz_b)
    add_arrows_xx(q, GG, min_edge_frac=min_edge_frac, r=r, arw_lng=arw_lng, lwd=lwd,
                  lamp=lamp, arw_clr=arw_clr, seg_clr=seg_clr)
    add_points_qq(q, GG, max_pt_sz = mx_pt_sz_q)
  })
  return(invisible())
}
