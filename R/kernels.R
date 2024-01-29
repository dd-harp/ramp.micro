
#' Make an exponential function to weight points by distance
#'
#' @param k decay by distance
#' @param s a scale parameter
#' @param gamma a shape parameter
#'
#' @return a function
#' @export
make_kF_exp = function(k=1, s=2, gamma=1){
  return(function(dd, w=1){
    wij = w*(exp(-k*(dd/s)^gamma))
    wij/max(wij)
  })
}

#' Make an power function to weight points by distance
#'
#' @param delta the power on distance
#' @param s a shape function
#'
#' @return a function
#' @export
make_kF_pwr = function(delta=1, s=1){
  return(function(dd, w=1){
    wij = w/(dd+s)^delta
    wij/max(wij)
  })
}

#' Make an power function that combines an exponential and power function
#'
#' @param p the weight on the power function
#' @param k decay by distance
#' @param s a scale parameter
#' @param gamma a shape parameter
#' @param delta the power on distance
#'
#' @return a function
#' @export
make_kF_mix = function(p=0.001, k=1, s=2, gamma=1, delta=1){
  return(function(dd, w=1){
    wij = (1-p)*w*(exp(-k*(dd/s)^gamma)) + p*w/(dd+s)^delta
    wij/max(wij)
  })
}
