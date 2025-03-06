
#' Draw arrows from one point set to another
#'
#' @param xy_launch origin point set
#' @param xy_land destination point set
#' @param M the dispersal matrix
#' @param min_edge_frac the fraction of mass to plot
#' @param r the radius of a ring around destination points
#' @param arw_lng the arrow length
#' @param lwd scale the line width
#' @param lamp arrow width scaling factor
#' @param clr the arrow color
#'
#' @return invisible(NULL)
#' @export
add_arrows_xy = function(xy_launch, xy_land, M, min_edge_frac=0.95,
                         r=0, arw_lng=0.1, lwd=2, lamp=1, clr="darkgreen"){
  n1 = dim(xy_launch)[1]
  n2 = dim(xy_land)[1]
  ij = edgeSubset(M, min_edge_frac)
  Mwt = lwd*M/max(M)
  wt = Mwt[ij]
  n = length(wt)
  if(length(clr)>1) clr_n = clr[ij[,1]] else clr_n = clr
  sapply(c(1:n), add_nth_arrow_xy, ij=ij, xy_launch=xy_launch, xy_land=xy_land, p=1, r=r, wt=wt, lamp=lamp, clr=clr_n, arw_lng=arw_lng) -> dmb
  return(invisible())
}

#' Draw arrows from a point set to itself
#'
#' @param xy the point set
#' @param M the dispersal matrix
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
add_arrows_xx = function(xy, M, min_edge_frac=.99, r=0, arw_lng=.1, lwd=5, lamp=1,
                         arw_clr = "darkgreen", seg_clr= "#CCCCCC"){
  Mtot = M + t(M)
  Mfrac = M/Mtot
  diag(Mtot) <- 0
  ij = edgeSubset(Mtot, min_edge_frac)
  Mwt = lwd*Mtot/max(Mtot)
  p = Mfrac[ij]
  ix = which(p>0.5)
  ij = ij[ix,]
  p = p[ix]
  wt = Mwt[ij]
  sapply(1:length(p), add_nth_segment, ij=ij, xy=xy, wt=wt, lamp=lamp, clr=seg_clr) -> dmb
  sapply(1:length(p), add_nth_arrow_xx, ij=ij, xy=xy, p=p, r=r, arw_lng=arw_lng, wt=wt, lamp=lamp, clr=arw_clr)-> dd
  return(invisible())
}

#' Add the nth arrow from a list of locations and indices
#'
#' @param n the arrow to plot
#' @param ij the indices of the arrows to plot
#' @param xy_launch the xy location of the search origin
#' @param xy_land the xy location of the search end
#' @param p fraction of length to plot
#' @param r the radius of a ring around destination points
#' @param arw_lng the arrow length
#' @param wt the line weight
#' @param lamp arrow width scaling factor
#' @param clr the arrow color
#'
#' @return invisible(NULL)
#' @export
add_nth_arrow_xy = function(n, ij, xy_launch, xy_land, p=1, r=0, arw_lng=0.1, wt=1, lamp=1, clr="#CCCCCC"){
  j = ij[n,1]; i=ij[n,2]
  x_i = xy_launch[i,1]; y_i=xy_launch[i,2]
  x_j = xy_land[j,1]; y_j=xy_land[j,2]
  if(length(clr)==1) clr_n=clr else clr_n=clr[n]
  if(length(p)==1) p_n=p else p_n=p[n]
  if(length(wt)==1) wt_n=wt else wt_n=wt[n]
  add_one_arrow(x_i, y_i, x_j, y_j, p_n, r, arw_lng, wt_n, lamp, clr_n)
  return(invisible())
}

#' Plot the n^th arrow
#'
#' @param n the arrow to plot
#' @param ij the indices of the arrows to plot
#' @param xy the xy locations
#' @param p fraction of length to plot
#' @param r the radius of a ring around destination points
#' @param arw_lng the arrow length
#' @param wt the line weight
#' @param lamp arrow width scaling factor
#' @param clr the arrow color
#'
#' @return invisible(NULL)
#' @export
add_nth_arrow_xx = function(n, ij, xy, p=1, r=0, arw_lng=0.1, wt=1, lamp=1, clr="#CCCCCC"){
  j = ij[n,1]; i=ij[n,2]
  x_i = xy[i,1]; y_i=xy[i,2]
  x_j = xy[j,1]; y_j=xy[j,2]
  if(length(clr)==1) clr_n=clr else clr_n=clr[n]
  if(length(p)==1) p_n=p else p_n=p[n]
  if(length(wt)==1) wt_n=wt else wt_n=wt[n]
  add_one_arrow(x_i, y_i, x_j, y_j, p_n, r, arw_lng, wt_n, lamp, clr_n)
  return(invisible())
}

#' Draw an arrow from (x_i, y_i) to (x_j, y_j)
#'
#' @param x_i source x value
#' @param y_i source y value
#' @param x_j destination x value
#' @param y_j destination y value
#' @param p fraction of length to plot
#' @param r the radius of a ring around destination points
#' @param arw_lng the arrow length
#' @param wt the line weight
#' @param lamp arrow width scaling factor
#' @param clr the arrow color
#'
#' @return invisible(NULL)
#' @export
add_one_arrow = function(x_i, y_i, x_j, y_j, p=1, r=0, arw_lng=0.1, wt=1, lamp=1, clr="#CCCCCC"){
  x0 = x_j + p*(x_i-x_j)
  y0 = y_j + p*(y_i-y_j)
  m = atan((y_i-y_j)/(x_i-x_j))
  m1=ifelse(length(m)>1, m[1], m)
  m2=ifelse(length(m)>1, m[2], m)
  x1 = x_j + sign(x_i-x_j)*r*cos(m)
  y1 = y_j + sign(x_i-x_j)*r*sin(m)
  if(sqrt((x0-x1)^2 + (y0-y1)^2) > arw_lng)
    suppressWarnings(graphics::arrows(x0, y0, x1, y1, length=arw_lng, lwd=lamp*log(1+wt), col=clr, lend=2))
  return(invisible())
}

#' Add the nth segment
#'
#' @param n an index for
#' @param ij the indices of the arrows to plot
#' @param xy the xy locations
#' @param wt the line weight
#' @param lamp arrow width scaling factor
#' @param clr the arrow color
#'
#' @return invisible(NULL)
#' @export
add_nth_segment = function(n, ij, xy, wt=1, lamp=1, clr="#CCCCCC"){
  i = ij[n,1]; j=ij[n,2]
  x_i = xy[i,1]; y_i=xy[i,2]
  x_j = xy[j,1]; y_j=xy[j,2]
  if(length(clr)==1) clr_n=clr else clr_n=clr[n]
  if(length(wt)==1) wt_n=wt else wt_n=wt[n]
  graphics::segments(x_i, y_i, x_j, y_j, lwd=lamp*log(1+wt_n), col=clr_n, lend="square", ljoin="bevel")
  return(invisible())
}
