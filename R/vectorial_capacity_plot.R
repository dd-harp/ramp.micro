
#' Visualize lifetime transmission from a site, per mosquito
#'
#' @param model a compound [list] defining a model
#' @param max_pt_sz_b the maximum cex for blood feeding sites
#' @param max_pt_sz_q the maximum cex for egg laying sites
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
plot_dispersal_V = function(model,
                    max_pt_sz_b = 0.7, max_pt_sz_q = 2,
                    min_edge_frac = 0.01, r=.01, arw_lng=0.05, lwd=2, lamp=1,
                    arw_clr="darkolivegreen4", seg_clr="orangered3"){
  with(model,{
    par(mar=c(2,2,2,2))
    frame_bq(b, q, mtl = "Potential Parasite Dispersal, per Mosquito")
    add_points_q(q, max_pt_sz = max_pt_sz_q)
    add_arrows_xx(b, V, min_edge_frac=min_edge_frac, r=r, arw_lng=arw_lng, lwd=lwd,
                  lamp=lamp, arw_clr=arw_clr, seg_clr=seg_clr)
    add_points_bb(b, V, max_pt_sz = max_pt_sz_b)
  }
  )}

#' Visualize lifetime transmission from a site by mosquito populations
#'
#' @param model a compound [list] defining a model
#' @param max_pt_sz_b the maximum cex for blood feeding sites
#' @param max_pt_sz_q the maximum cex for egg laying sites
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
plot_dispersal_VV = function(model,
                     max_pt_sz_b = 0.7, max_pt_sz_q = 2,
                     min_edge_frac = 0.01, r=.01, arw_lng=0.05, lwd=2, lamp=1,
                     arw_clr="springgreen4", seg_clr="firebrick3"){
  with(model,{
    par(mar=c(2,2,2,2))
    frame_bq(b, q, mtl = "Potential Parasite Dispersal, Population")
    add_points_q(q, max_pt_sz = max_pt_sz_q)
    add_arrows_xx(b, VV, min_edge_frac=min_edge_frac, r=r, arw_lng=arw_lng, lwd=lwd,
                  lamp=lamp, arw_clr=arw_clr, seg_clr=seg_clr)
    add_points_bb(b, VV, max_pt_sz = max_pt_sz_b)
  })
  return(invisible())
}
