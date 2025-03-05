
#' Make an exponential function for weight by distance
#'
#' @description This returns a function of the form
#' \deqn{F_w (d, \omega=1) = \omega_j e^{-k \left( \frac{d_{i,j}}{s}\right)^\gamma}}
#' where \eqn{s} and \eqn{\gamma} are shape parameters, \eqn{k} is the rate
#' parameter, and \eqn{\omega} is a weight.
#'
#' In effect, \eqn{s} is the location of a shoulder, and for \eqn{\gamma>1}, the decay is
#' slower for \eqn{d<s}.
#'
#' The function returned accepts \eqn{\omega} as
#' an optional argument so that it can be passed at the time of simulation.
#'
#' By default, the function returns scaled values -- the maximum is 1.
#'
#' @param k decay by distance
#' @param s a scale parameter
#' @param gamma a shape parameter
#'
#' @return a function
#' @export
#' @examples
#' kF1 = make_kF_exp(k=1, s=1, gamma=1.5)
#' kF2 = make_kF_exp(k=2, s=0.1, gamma=2)
#' dd = seq(0, 2, by = 0.01)
#' plot(dd, kF1(dd), type = "l", ylab = "Weight", xlab = "Distance")
#' lines(dd, kF2(dd))
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
#' @param s1 a scale parameter
#' @param s2 a scale parameter
#' @param gamma a shape parameter
#' @param delta the power on distance
#'
#' @return a function
#' @export
make_kF_mix = function(p=0.001, k=1, s1=1, s2=1, gamma=1, delta=1){
  return(function(dd, w=1){
    wij = (1-p)*w*(exp(-k*(dd/s1)^gamma)) + p*w/(dd+s2)^delta
    wij/max(wij)
  })
}
