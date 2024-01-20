


plotGraph = function(mod, net, cut=NULL,  alg = "wt", clr=grey(0.7), adj=2, pw=1, mtl = "", stretch=0.1, lwd=2){with(mod,{
  if(alg == "wt") clusters = net$clusters_walktrap
  if(alg == "gr") clusters = net$clusters_greedy
  memix = if(is.null(cut)){
    membership(clusters)
  } else {
    cut_at(clusters, cut)
  }
  nC = max(memix)
  clrs = turbo(nC)
  frame_bq(b,q,mtl)
  if(net$type == "b"){
    xy = b
    arrowsX2X(b, net$M, clr = clr, adj=adj)
    addP.b(b, rowSums(net$M), clr=clrs[memix], adj=adj, pw=pw)
  }
  if(net$type == "q"){
    xy = q
    arrowsX2X(q, net$M, clr = clr, adj=adj)
    addP.q(q, rowSums(net$M), clr=clrs[memix], adj=adj, pw=pw)
  }
  if(net$type == "q"){
    xy = rbind(b,q)
    arrowsX2X(xy, net$M, clr = clr, adj=adj)
    wts = rowSums(net$M)
    addP.b(b, wts[1:nb], clr=clrs[memix[1:nb]], adj=adj, pw=pw)
    addP.q(q, wts[nb+1:nq], clr=clrs[memix[nb+1:nq]], adj=adj, pw=pw)
  }
  addConvexHulls(memix, xy, clrs, stretch, lwd)
})}

tilesGraph = function(mod, net, cut=NULL,  alg = "wt", clrs=NULL, adj=2, pw=1, mtl = "", stretch=0.1, lwd=2){with(mod,{
  if(alg == "wt") clusters = net$clusters_walktrap
  if(alg == "gr") clusters = net$clusters_greedy
  memix = if(is.null(cut)){
    membership(clusters)
  } else {
    cut_at(clusters, cut)
  }
  nC = max(memix)
  clrs = turbo(nC)
  if(is.null(clrs)) clrs = viridis(nC)
  if(net$type == "b"){
    plot(tilesB, fillcol = clrs[memix])
  }
  if(net$type == "q"){
    plot(tilesQ, fillcol = clrs[memix])
  }
})}

communityArrows = function(mod, i, cut=NULL, clrP=turbo, mnwd=0.05, bbend=1, adj=2, endd=0.75){
  net = getNet.i(mod, i)
  clusters = net$clusters_walktrap
  M = getM.i(mod, i)
  if(is.null(cut)){
    memix = membership(clusters)
  } else{
    memix = cut_at(clusters, cut)
  }
  pal = clrP(max(memix))
  clrs = pal[memix]
  if(net$type == "b") xy = mod$b
  if(net$type == "q") xy = mod$q
  with(mod, frame_bq(b, q))
  #  arrowsX2X(xy,M, mnwd, bbend, adj, endd, clr=clrs)
  arrowsX2X_fm(xy, M, .9, bbend, adj, endd, clr=clrs)
}

