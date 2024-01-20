
getConvexHull.i = function(i, memix, xy){
  ixj = which(memix == i)
  hpts <- chull(xy[ixj,])
  ixk = c(hpts, hpts[1])
  return(xy[ixj[ixk],])
}

modConvexHulls.i = function(mod, i, cut=NULL, clr_scheme=turbo, stretch=0.1, lwd=2){
  net = getNet.i(mod, i)
  modConvexHulls(mod, net, cut, clr_scheme, stretch, lwd)
}

modConvexHulls = function(mod, net, cut=NULL, clrP=turbo, stretch=0.1, lwd=2, clrs=NULL){
  clusters = net$clusters_walktrap
  memix = if(is.null(cut)){
    membership(clusters)
  } else {
    cut_at(clusters, cut)
  }
  if(net$type == "b"){
    xy = mod$b
  }
  if(net$type == "q"){
    xy = mod$q
  }
  if(net$type == "bq"){
    xy = with(mod,rbind(b, q))
  }
  if(is.null(clrs)) clrs = clrP(max(memix))
  addConvexHulls(memix, xy, clrs, stretch, lwd)
}

addConvexHulls = function(memix, xy, clrs, stretch=0.1, lwd=2){
  for(i in 1:max(memix)){
    hxy = getConvexHull.i(i, memix, xy)
    sxy = stretchHull(hxy, 1+stretch)
    polygon(sxy[,1], sxy[,2], border=clrs[i], lwd=lwd)
  }
}

stretchHull = function(xy, fac){
  cx = mean(xy[,1])
  cy = mean(xy[,2])
  xn = xy[,1] - cx
  yn = xy[,2] - cy

  #x = r*cos(theta)
  #y = r*sin(theta)
  r = sqrt(xn^2 + yn^2)
  theta = atan2(yn,xn)
  sxy =  cbind(
    r*fac*cos(theta) + cx,
    r*fac*sin(theta) + cy
  )
  return(sxy)
}


