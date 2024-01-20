
plotMeta = function(mod, i, cut=NULL, lwd=2, bbend=3, mtl = NULL){
  net = getNet.i(mod, i)
  meta = net2meta(mod, i, cut)
  n = dim(meta$centers)[1]
  clrs = turbo(n)[sample(1:n)]
  with(mod, frame_bq(b,q, mtl))
  points(meta$centers, col = clrs)
  modConvexHulls(mod, net, cut, lwd=lwd, clrs=clrs)
  bentArrowsX2X(meta$centers, meta$p.connect, bbend=bbend, clr=clrs)
}



triplePlotHulls = function(mod, cut=NULL, pairs="Kxy", fm = 0.9){with(mod,{
  if(pairs == "Kxy"){
    Mb = Kbb
    Mq = Kqq
    bnet = Kbb_net
    qnet = Kqq_net
    mtlb = expression(K[b%->%b])
    mtlq = expression(K[q%->%q])
  }
  if(pairs == "VG"){
    Mb = VC
    Mq = GG
    bnet = VC_net
    qnet = GG_net
    mtlb = expression("Potential Parasite Dispersal, Population")
    mtlq = expression("Lifetime Egg Dispersal, Population")
  }
  if(is.null(cut)){
    membershipIXb = membership(bnet$clusters_walktrap)
    Nb = max(membershipIXb)
    membershipIXq = membership(qnet$clusters_walktrap)
    Nq = max(membershipIXq)
  }
  else{
    Nb=Nq=cut
    membershipIXb = cut_at(bnet$clusters_walktrap, N)
    membershipIXq = cut_at(qnet$clusters_walktrap, N)
  }
  palb = turbo(Nb)[sample(1:Nb)]
  clrsb = palb[membershipIXb]
  palq = turbo(Nq)[sample(1:Nq)]
  clrsq = palq[membershipIXq]

  frame_bq(b, q, mtl = mtlb)
  arrowsX2Y_fm(b, b, Mb, clr=clrsb, bbend=2, fracMass=fm)
  addP.b(b, wts = rowSums(Mb), adj=3, clr=clrsb)

  frame_bq(b, q, mtl = expression("Convex Hulls"))
  addConvexHulls(membershipIXb, b, palb, stretch=-0.1, lwd=4)
  addConvexHulls(membershipIXq, q, palq, stretch=-0.1, lwd=4)

  frame_bq(b, q, mtl = mtlq)
  arrowsX2Y_fm(q, q, Mq, clr=clrsq, bbend=2, fracMass=fm)
  addP.q(q, wts = rowSums(Mq), adj=3, clr=clrsq)
})}
