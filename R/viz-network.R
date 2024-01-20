



cutoffValue = function(M, fracMass=0.9){
  # returns a cutoffValue value
  vals = sort(M, decreasing=TRUE)
  max(which(cumsum(vals)<fracMass*sum(M))) -> mx
  return(vals[mx])
}

edgeSubset = function(M, fracMass=0.9){
  # returns the edges that make up
  # fracMass of all flows
  V = cutoffValue(M, fracMass)
  which(M >= V, arr.ind=T)
}

flow = function(xi, yi, xj, yj, p=1, r=0, wt=1, clr=grey(0.7), lng=0.25, lamp=1){
  # draw an arrow from (xi, yi) to (xj, yj)
  x1 = xj + p*(xi-xj)
  y1 = yj + p*(yi-yj)
  m = atan((yj-yi)/(xj-xi))
  x2 = xj - r*cos(m)
  y2 = yj - r*sin(m)
  arrows(x1, y1, x2, y2, length=lng, lwd=lamp*log(1+wt), col=clr, lend=2)
  # segments(xi, yi, xj, yj, lwd=wt, col=clr, lend=2)
}

# to plot flows
flow.i = function(n, ij, xy, p, r, wt, clr, lng){
  i = ij[n,1]; j=ij[n,2]
  xi = xy[i,1]; yi=xy[i,2]
  xj = xy[j,1]; yj=xy[j,2]
  if(length(clr)==1) clr else clr[n]
  flow(xj,yj,xi,yi,p[n],r,wt[n],clr,lng)
}

flow2.i = function(n, ij, xy1, xy2, wt, clr, lng){
  j = ij[n,1]; i=ij[n,2]
  xi = xy1[i,1]; yi=xy1[i,2]
  xj = xy2[j,1]; yj=xy2[j,2]
  if(length(clr)==1) clr0=clr else clr0=clr[n]
  flow(xi, yi, xj, yj, wt=wt[n], clr=clr0, lng=lng)
  # segments(xi, yi, xj, yj, lwd=wt[n], col=clr, lend=2)
}

edge = function(xi, yi, xj, yj, wt, clr){
  # draw a segment between (xi, yi) and (xj, yj)
  segments(xi, yi, xj, yj, lwd=wt, col=clr, ljoin="bevel")
}

edge.i = function(n, ij, xy, wt, clr, lamp=1){
  i = ij[n,1]; j=ij[n,2]
  xi = xy[i,1]; yi=xy[i,2]
  xj = xy[j,1]; yj=xy[j,2]
  segments(xi, yi, xj, yj, lwd=lamp*log(1+wt[n]), col=clr, lend="square", ljoin="bevel")
}



arrowsX2Y = function(xy1, xy2, M, cutat=0.99, lng=0.02, adj=2, clr="red"){
  n1 = dim(xy1)[1]
  n2 = dim(xy2)[1]
  ij = edgeSubset(M, cutat)
  Mwt = adj*M/max(M)
  wt = Mwt[ij]
  n = length(wt)
  if(length(clr)>1) clr = clr[ij[,1]]
  sapply(c(1:n), flow2.i, ij=ij, xy1=xy1, xy2=xy2, wt=wt, clr=clr, lng=lng) -> dmb
}

bentArrowsX2Y = function(xy1, xy2, M, mnwd=0.05, bbend=0, endd=0.75, adj=2, clr="red"){
  n1 = dim(xy1)[1]
  n2 = dim(xy2)[1]
  if(length(clr == 1)) clr = rep(clr, n1)
  adj = adj/max(M)
  for(i in 1:n1)
    for (j in c(1:n2))
      if (M[j,i] > mnwd*max(M)){
        fac = M[j,i]*adj
        curvedarrow(xy1[i,], xy2[j,],
                    segment=c(0.1, endd), lwd=fac, curve = 0.01*bbend,
                    arr.pos=endd, arr.length=0.15*fac, arr.width=0.1*fac,
                    arr.type="curved", arr.col = clr, lcol = clr[i])
      }
}



arrowsX2Xv1 = function(xy, M, cutat=.99, r=0,
                       clrA="red", clrS = grey(0.6),
                       fac1=5, fac2=1, lng=.1, wd=.015){
  Mtot = M + t(M)
  Mfrac = M/Mtot
  diag(Mtot) <- 0
  ij = edgeSubset(Mtot, cutat)
  Mwt = fac1*Mtot/max(Mtot)
  p = Mfrac[ij]
  ix = which(p>0.5)
  ij = ij[ix,]
  p = p[ix]
  wt = Mwt[ij]
  sapply(1:length(p), edge.i, ij=ij, xy=xy, wt=wt, clr=clrS) -> dmb
  sapply(1:length(p), flow.i, ij=ij, xy=xy, wt=wt, p=p, r=r, clr=clrA, lng=lng)-> dd
}

arrowsX2Xv2 = function(xy, M, cutat=.99, r=0,
                       clrA="red", clrS = grey(0.6),
                       fac1=5, fac2=1, lng=.05, wd=.015){
  Mtot = M + t(M)
  Mfrac = M/Mtot
  ij = edgeSubset(Mtot, cutat)
  ix = which(ij[,1] == ij[,2])
  ij = ij[-ix,]
  Mwt = fac1*Mtot/max(Mtot)
  p = Mfrac[ij]
  ix = which(p>0.5)
  ij = ij[ix,]
  p = p[ix]
  wt = Mwt[ij]
  sapply(1:length(p), edge.i, ij=ij, xy=xy, wt=wt, clr=clrS) -> dmb
  sapply(1:length(p), flow.i, ij=ij, xy=xy, wt=wt, p=p, r=r, clr=clrA, lng=lng)-> dd
}

arrowsX2X = arrowsX2Xv1

bentArrowsX2X = function(xy, M, mnwd=0.05, bbend=1, adj=2, endd=0.75, clr = "red" ){
  n = dim(xy)[1]
  if(length(clr == 1)) clr = rep(clr, n)
  diag(M) <- 0
  adj = adj/max(M)
  for(i in 1:n)
    for(j in c(1:n)[-i])
      if (M[j,i]> mnwd*max(M)){
        fac = M[j,i]*adj
        curvedarrow(xy[i,], xy[j,],
                    segment=c(0.1, endd), lwd=fac, curve = 0.01*bbend,
                    arr.pos=endd, arr.length=0.15*fac, arr.width=0.1*fac,
                    arr.type="curved", arr.col = clr, lcol = clr[i])
      }
}

decompM=function(M){
  diagonal = diag(M)
  sym = pmin(M, t(M))
  flow = M-sym
  list(sym=sym, flow=flow, diagonal=diagonal)
}
