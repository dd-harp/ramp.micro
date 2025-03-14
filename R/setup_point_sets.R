#' Make a random point set, xy
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

#' Set up a lattice with n x n points
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

#' Use a set of `xy` values as seeds for clusters
#'
#' @param xy a set of seeds
#' @param nc the number of points per cluster
#' @param vr the variance in the dispersion around the seed
#'
#' @return a set of xy locations
#' @export
clusters_xy = function(xy, nc=1, vr=1){
  np = dim(xy)[1]
  nc = stats::rpois(np, nc)
  for(i in 1:np){
    if(nc[i]>0){
      xi = xy[i,1] + stats::rnorm(nc[i], 0, vr)
      yi = xy[i,2] + stats::rnorm(nc[i], 0, vr)
      xyi = cbind(xi,yi)
      xy = rbind(xy, xyi)
    }
  }
  xy
}

#' Setup n clusters with m points per cluster
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
    clust = 1+stats::rpois(1,m-1)
    xi = xy[i,1] + stats::rnorm(clust, 0, vr)
    yi = xy[i,2] + stats::rnorm(clust, 0, vr)
    xyi = cbind(xi,yi)
    xy = rbind(xy, xyi)
  }
  xy
}
