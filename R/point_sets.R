#' Make a random xy point set
#'
#' @param n the number of points
#' @param mn the lower limit for location values
#' @param mx the upper limit for location values
#'
#' @return a set of xy locations
#' @export
unif_xy = function(n, mn, mx){
  x = stats::runif(n, mn, mx)
  y = stats::runif(n, mn, mx)
  cbind(x,y)
}

#' Title
#'
#' @param n the number of points
#' @param mn the lower limit for location values
#' @param mx the upper limit for location values
#'
#' @return a set of xy locations
#' @export
lattice = function(n, mn, mx){
  points = seq(mn,mx,length.out=n)
  m = matrix(points, n, n)
  cbind(as.vector(m), as.vector(t(m)))
}

#' Use the xy values as the seed locations to generate clusters
#'
#' @param xy a set of seeds
#' @param nc the number of points per cluster
#' @param vr the variance in the dispersion around the seed
#'
#' @return a set of xy locations
#' @export
clusters_xy = function(xy, nc=1, vr=1){
  np = dim(xy)[1]
  if(length(nc) < np)
    nc = 1 + stats::rpois(nc, 1)
  for(i in 1:np){
    xi = xy[i,1] + stats::rnorm(nc[i], 0, vr)
    yi = xy[i,2] + stats::rnorm(nc[i], 0, vr)
    xyi = cbind(xi,yi)
    xy = rbind(xy, xyi)
  }
  xy
}

#' Title
#'
#' @param n the number of seed points points
#' @param m the number of points per cluster
#' @param vr the variance in the dispersion around the seed
#' @param mn the lower limit for location values
#' @param mx the upper limit for location values
#'
#' @return a set of xy locations
#' @export
clusters_nm = function(n, m, vr, mn, mx){
  xy = unif_xy(n, mn, mx)
  for(i in 1:n){
    clust = 1+rpois(1,m-1)
    xi = xy[i,1] + stats::rnorm(clust, 0, vr)
    yi = xy[i,2] + stats::rnorm(clust, 0, vr)
    xyi = cbind(xi,yi)
    xy = rbind(xy, xyi)
  }
  xy
}
