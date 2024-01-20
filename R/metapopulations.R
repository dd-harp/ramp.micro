

getM.i = function(mod, i){
  if(i==1) M = mod$Kbb
  if(i==2) M = mod$Kqq
  if(i==3) M = mod$G
  if(i==4) M = mod$GG
  if(i==5) M = mod$V
  if(i==6) M = mod$VC
  if(i==7) M = mod$bigM
  if(i==8) M = mod$bigMM
  return(M)
}

getNet.i = function(mod, i){
  if(i==1) net = mod$Kbb_net
  if(i==2) net = mod$Kqq_net
  if(i==3) net = mod$G_net
  if(i==4) net = mod$GG_net
  if(i==5) net = mod$V_net
  if(i==6) net = mod$VC_net
  if(i==7) net = mod$M_net
  if(i==8) net = mod$MM_net
  return(net)
}

net2meta =function(mod, i, cut=NULL){
  net = getNet.i(mod, i)
  clusters = net$clusters_walktrap
  if(net$type == "b") xy = mod$b
  if(net$type == "q") xy = mod$q
  M = getM.i(mod, i)
  if(is.null(cut)){
    memix = membership(clusters)
  } else{
    memix = cut_at(clusters, cut)
  }
  dm = max(memix)
  pop = matrix(0, dm, dm)
  for(i in 1:dm){
    for(j in 1:dm){
      ixi = which(memix==i)
      ixj = which(memix==j)
      pop[i,j] = sum(M[ixi,ixj])
    }
  }
  cntr = c()
  for( i in 1:dm){
    ix = which(memix==i)
    cx = mean(xy[ix,1])
    cy = mean(xy[ix,2])
    cxy = cbind(cx,cy)
    cntr = rbind(cntr, cxy)
  }
  self = diag(pop)
  p.self = self/colSums(pop)
  pop1 = pop
  diag(pop1)<-0
  p.connect = pop1%*% diag(1/colSums(pop))
  list(pop=pop, self=self, p.self=p.self, p.connect=p.connect, connect=pop1, centers=cntr)
}
