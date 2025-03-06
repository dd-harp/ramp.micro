
#' Draw bent arrows from one point set to another
#'
#' @param xy_launch origin point set
#' @param xy_land destination point set
#' @param M the dispersal matrix
#' @param mnwd the minimum width to plot
#' @param bbend a parameter to bend arrows
#' @param endd relative position of arrowhead
#' @param adj set the maximum cex for points
#' @param clr the color
#'
#' @return invisible(NULL)
#' @export
add_bent_arrows_xy = function(xy_launch, xy_land, M, mnwd=0.05, bbend=0, endd=0.75, adj=2, clr="red"){
  n1 = dim(xy_launch)[1]
  n2 = dim(xy_land)[1]
  if(length(clr == 1)) clr = rep(clr, n1)
  adj = adj/max(M)
  for(i in 1:n1)
    for (j in c(1:n2))
      if (M[j,i] > mnwd*max(M)){
        fac = M[j,i]*adj
        diagram::curvedarrow(xy_launch[i,], xy_land[j,],
                    segment=c(0.1, endd), lwd=fac, curve = 0.01*bbend,
                    arr.pos=endd, arr.length=0.15*fac, arr.width=0.1*fac,
                    arr.type="curved", arr.col = clr, lcol = clr[i])
      }
  return(invisible())
}


#' Title
#'
#' @param xy the point set
#' @param M the dispersal matrix
#' @param mnwd the minimum width to plot
#' @param bbend a parameter to bend arrows
#' @param endd relative position of arrowhead
#' @param adj set the maximum cex for points
#' @param clr the color
#'
#' @return invisible(NULL)
#' @export
add_bent_arrows_xx = function(xy, M, mnwd=0.05, bbend=1, adj=2, endd=0.75, clr = "red"){
  n = dim(xy)[1]
  if(length(clr == 1)) clr = rep(clr, n)
  diag(M) <- 0
  adj = adj/max(M)
  for(i in 1:n)
    for(j in c(1:n)[-i])
      if (M[j,i]> mnwd*max(M)){
        fac = M[j,i]*adj
        diagram::curvedarrow(xy[i,], xy[j,],
                             segment=c(0.1, endd), lwd=fac, curve = 0.01*bbend,
                             arr.pos=endd, arr.length=0.15*fac, arr.width=0.1*fac,
                             arr.type="curved", arr.col = clr, lcol = clr[i])
      }
  return(invisible())
}
