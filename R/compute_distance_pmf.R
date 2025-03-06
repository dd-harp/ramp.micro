
#' Compute distances dispersed
#'
#' @param xy locations of a set of points
#' @param M a dispersal matrix
#'
#' @return a [list] with the PMF, CMF, mean distance dispersed, and the distance matrix
#' @export
dispersal_PMF = function(xy, M){
  dxy = as.matrix(stats::dist(xy, diag=TRUE, upper=TRUE))
  dxy = as.vector(dxy)
  M = as.vector(M)
  ot = order(dxy)
  dxy = dxy[ot]
  PMF = M[ot]/sum(M)
  CMF = cumsum(PMF)
  list(mean = sum(dxy*PMF), pmf=PMF, cmf=CMF, dxy=dxy)
}

#' Plot the PMF for distances dispersed
#'
#' @param xMFobj a [list] returned from [dispersalPMF]
#' @param mtl a title for the plot
#' @param clr the colors
#' @param thresh a threshold for
#'
#' @return invisible(NULL)
#' @export
plotDDpmf = function(xMFobj, mtl=NULL, clr="black", thresh=.99){
  with(xMFobj,{
    ix = 1:min(which(cmf>thresh))
    dxy = dxy[ix]
    pmf = pmf[ix]
    plot(dxy, pmf, main = mtl, type ="h", xlab = "Distance", ylab = "PMF", lwd=2)
    segments(mean, 0, mean, max(pmf), col = clr, lwd=2, lty=2)
  })
  return(invisible())
}

#' Plot the CMF for distances dispersed
#'
#' @param xMFobj a [list] returned from [dispersalPMF]
#' @param mtl a title for the plot
#' @param clr the colors
#' @param thresh a threshold for
#'
#' @return invisible(NULL)
#' @export
plotDDcmf = function(xMFobj, mtl=NULL, clr="black", thresh=.99){
  with(xMFobj,{
    ix = 1:min(which(cmf>thresh))
    dxy = dxy[ix]
    cmf = cmf[ix]
    plot(dxy, cmf, main = mtl, type ="l", xlab = "Distance", ylab = "CMF", ylim = c(0,1), lwd=2)
    segments(mean, 0, mean, 1, col = clr, lwd=2, lty=2)
  })
  return(invisible())
}
