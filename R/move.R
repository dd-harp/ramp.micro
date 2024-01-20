
dispersalPMF = function(xy, M){
  dxy = as.matrix(dist(xy, diag=TRUE, upper=TRUE))
  dxy = as.vector(dxy)
  M = as.vector(M)
  ot = order(dxy)
  dxy = dxy[ot]
  PMF = M[ot]/sum(M)
  CMF = cumsum(PMF)
  list(dxy=dxy, pmf=PMF, cmf=CMF, mean = sum(dxy*PMF))
}

plotDDpmf = function(xMFobj, mtl=NULL, clr="black", thresh=.99){with(xMFobj,{
  ix = 1:min(which(cmf>thresh))
  dxy = dxy[ix]
  pmf = pmf[ix]
  plot(dxy, pmf, main = mtl, type ="h", xlab = "Distance", ylab = "PMF", lwd=2)
  segments(mean, 0, mean, max(pmf), col = clr, lwd=2, lty=2)
})}

plotDDcmf = function(xMFobj, mtl=NULL, clr="black", thresh=.99){with(xMFobj,{
  ix = 1:min(which(cmf>thresh))
  dxy = dxy[ix]
  cmf = cmf[ix]
  plot(dxy, cmf, main = mtl, type ="l", xlab = "Distance", ylab = "CMF", ylim = c(0,1), lwd=2)
  segments(mean, 0, mean, 1, col = clr, lwd=2, lty=2)
})}
