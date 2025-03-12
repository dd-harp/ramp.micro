# Visualize dispersal among point sets

#' Visualize a dispersal matrix among two point sets
#'
#' @param xy_launch origin point set
#' @param xy_land destination point set
#' @param M the dispersal matrix
#' @param mx_pt_sz1 the maximum cex for the launch point sest
#' @param mx_pt_sz2 the maximum cex for the landing point set
#' @param pt_clr1 point color for launch
#' @param pt_clr2 point color for landing
#' @param min_edge_frac the minimum fraction to plot edge
#' @param r the radius of a ring around destination points
#' @param arw_lng the arrow length
#' @param arw_lwd scale the line width
#' @param lamp arrow width scaling factor
#' @param arw_clr the arrow color
#' @param mtl = plot title
#'
#' @return no visible return value
#' @export
plot_matrix_xy = function(xy_launch, xy_land, M,
                    mx_pt_sz1=0.7, mx_pt_sz2=1.5, pt_clr1 = "salmon", pt_clr2="lightblue",
                    min_edge_frac=0.01, r=0, arw_lng=0.02, arw_lwd=2, lamp=1, arw_clr="tomato",
                    mtl = ""){
  plot(rbind(xy_launch, xy_land), type = "n", xaxt = "n", yaxt = "n",
       xlab = "", ylab = "", main = mtl)
  graphics::points(xy_launch, col=pt_clr1, cex=mx_pt_sz1, pch=15)
  w = rowSums(M)
  wts = mx_pt_sz2*w/max(w)
  add_arrows_xy(xy_launch, xy_land, M, min_edge_frac, r, arw_lng, arw_lwd, lamp, arw_clr)
  graphics::points(xy_land, col = pt_clr2, cex=wts, pch=19)
  return(invisible())
}

#' Visualize a dispersal matrix from one point set to itself
#'
#' @param xy the point set
#' @param M the dispersal matrix
#' @param mx_pt_sz set cex for the largest point
#' @param pt_clr the point color
#' @param min_edge_frac the minimum fraction to plot edge
#' @param r the radius of a ring around destination points
#' @param arw_lng the arrow length
#' @param lwd set the maximum line width
#' @param lamp arrow width scaling factor
#' @param arw_clr the color to draw the arrow (asymetric part)
#' @param seg_clr the color to draw the segment (symmetric part)
#' @param mtl = plot title
#'
#' @return no visible return value
#' @export
plot_matrix_xx = function(xy, M,
                    mx_pt_sz=2, pt_clr = "tomato",
                    min_edge_frac=0.01, r=0.01, arw_lng=0.05, lwd=1, lamp=1,
                    arw_clr ="#e2739655", seg_clr = "#00000022",
                    mtl = ""){
  w = rowSums(M)
  wts = mx_pt_sz*w/max(w)
  plot(xy, xaxt = "n", yaxt = "n", col = pt_clr, pch = 19,
       xlab = "", ylab = "", main = mtl, cex = wts)
  add_arrows_xx(xy, M, min_edge_frac, r, arw_lng, lwd, lamp, arw_clr, seg_clr)
  graphics::points(xy, col = pt_clr, cex=wts, pch=19)
  return(invisible())
}



