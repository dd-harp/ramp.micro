

computeG= function(simObject, Tmax){
  UseMethod("computeG", simObject)
}

computeG.BQ = function(simObject, Tmax=50){with(simObject,{
  Q0 = diag(1, nq)
  B = Mqb %*% Q0
  Q = G = 0*Q0
  for(i in 1:Tmax){
    Bt = Mbb %*% B + Mqb %*% Q
    Qt = Mbq %*% B + Mqq %*% Q
    eggs = ova*psiQ*Q
    G = G + eggs
    B = Bt; Q=Qt
  }
  simObject$G = G
  return(simObject)
})}

computeG.BQS = function(simObject, Tmax=50){with(simObject,{
  Q0 = diag(1, nq)
  B = Mqb %*% Q0
  S = Mqs %*% Q0
  Q = G = 0*Q0
  for(i in 1:Tmax){
    Bt = Mbb %*% B + Mqb %*% Q + Msb%*%S
    Qt = Mbq %*% B + Mqq %*% Q + Msq%*%S
    St = Mbs %*% B + Mqs %*% Q + Mss%*%S
    eggs = ova*psiQ*Q
    G = G + eggs
    B = Bt; Q=Qt; S=St
  }
  simObject$G = G
  return(simObject)
})}

computeGG = function(simObject){with(simObject,{
  if(!exists("simObject$steadyState$Q")) simObject = steadyState(simObject)
  simObject$GG = with(simObject,G %*% diag(as.vector(steadyState$Q)))
  return(simObject)
})}

GProfile = function(
    mod, cutat=0.95, pw=0.5, adj=2, fac1=5, fac2=1,
    clrA="lightblue", clrS="salmon", clrQA="blue", clrQB="darkblue"
){with(mod,{
  par(mar=c(2,2,2,2))
  frame_bq(b, q, mtl = "Lifetime Egg Dispersal, per Mosquito")
  arrowsX2X(q, G, cutat=cutat, clrA=clrA, clrS=clrS, lng=0, fac1=fac1, fac2=fac2)
  addP.qq(q, G, pw=pw, adj=adj, colA=clrQA, colB=clrQB)
})}

GGProfile = function(
    mod, cutat=0.95, adj=2, pw=0.5, fac1=5, fac2=1,
    clrA="steelblue", clrS="chocolate", clrQA="blue", clrQB="darkblue"
){with(mod,{
  par(mar=c(2,2,2,2))
  frame_bq(b, q, mtl = "Lifetime Egg Dispersal, Population")
  arrowsX2X(q, GG, cutat=cutat, clrA=clrA, clrS=clrS, lng=0, fac1=fac1, fac2=fac2)
  addP.qq(q, GG, pw=pw, adj=adj, colA=clrQA, colB=clrQB)
})}

EGGProfile = function(mod, withHist=FALSE){with(mod,{
  if(withHist == TRUE) par(mfcol = c(2,2))
  if(withHist == FALSE) par(mfcol = c(1,2))

  frame_bq(b, q)
  arrowsX2X(q, G, colrA="blue")
  addP.qq(q, G, pw=0.5, adj=3)

  if(withHist == TRUE) hist(rowSums(G), col="blue")


  frame_bq(b, q)
  arrowsX2X(q, G, colrA="darkblue")
  addP.qq(q, GG, pw=0.5, adj=3)

  if(withHist == TRUE) hist(rowSums(GG), col="darkblue")
})}
