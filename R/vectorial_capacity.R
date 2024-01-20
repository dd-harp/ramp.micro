




computeV = function(simObject, Tmax){
  UseMethod("computeV", simObject)
}

computeV.BQ = function(simObject, Tmax=100){with(simObject,{
  Q = Mbq %*% diag(1, nb)
  B = diag(0, nb)

  for (i in 1:eip){
    Bt = Mbb %*% B + Mqb %*% Q
    Qt = Mbq %*% B + Mqq %*% Q
    B=Bt; Q=Qt
  }

  Vt = 0*B
  for (i in 1:Tmax){
    Vt = Vt + pB*psiB*B
    Bt = Mbb %*% B + Mqb %*% Q
    Qt = Mbq %*% B + Mqq %*% Q
    B=Bt; Q=Qt
  }
  simObject$V = Vt
  return(simObject)
})}

computeV.BQS = function(simObject, Tmax=100){with(simObject,{
  Q = Mbq %*% diag(1, nb)
  B = diag(0, nb)
  S = matrix(0, ns, nb)

  for (i in 1:eip){
    Bt = Mbb %*% B + Mqb %*% Q + Msb%*%S
    Qt = Mbq %*% B + Mqq %*% Q + Msq%*%S
    St = Mbs %*% B + Mqs %*% Q + Mss%*%S
    B=Bt; Q=Qt
  }

  Vt = 0*B
  for (i in 1:Tmax){
    Vt = Vt + pB*psiB*B
    Bt = Mbb %*% B + Mqb %*% Q + Msb%*%S
    Qt = Mbq %*% B + Mqq %*% Q + Msq%*%S
    St = Mbs %*% B + Mqs %*% Q + Mss%*%S
    B=Bt; Q=Qt
  }
  simObject$V = Vt
  return(simObject)
})}

computeVC = function(simObject, Tmax=100){with(simObject,{
  if(!exists("simObject$steadyState$B")) simObject = steadyState(simObject)
  simObject$VC = with(simObject, V %*% diag(as.vector(steadyState$B)))
  return(simObject)
})}

VProfile = function(
    mod, cutat=0.99, adj=2, pw=0.5, fac1=5, fac2=1,
    clrS="darkolivegreen4", clrA="orangered3",
    clrBA="darkolivegreen4", clrBB="orangered3"
){
  with(mod,{
    par(mar=c(2,2,2,2))
    frame_bq(b, q, mtl = "Potential Parasite Dispersal, per Mosquito")
    arrowsX2X(b, V, cutat=cutat, clrA=clrA, clrS=clrS, lng=0, fac1=fac1, fac2=fac2)
    addP.bb(b, V, pw=0.5, adj=adj, colA=clrBA, colB=clrBB)
  }
  )}

VCProfile = function(
    mod, cutat=0.99, adj=2, pw=0.5, fac1=5, fac2=1,
    clrS="springgreen4", clrA="firebrick3",
    clrBA="darkolivegreen4", clrBB="orangered3"
){
  with(mod,{
    par(mar=c(2,2,2,2))
    frame_bq(b, q, mtl = "Potential Parasite Dispersal, Population")
    arrowsX2X(b, VC, cutat=cutat, clrA=clrA, clrS=clrS, lng=0, fac1=fac1, fac2=fac2)
    addP.bb(b, VC, pw=pw, adj=adj, colA=clrBA, colB=clrBB)
  }
  )}

VVCProfile = function(mod, withHist=FALSE){
  if(withHist == TRUE) par(mfcol = c(2,2), mar = c(2,2,2,2))
  if(withHist == FALSE) par(mfcol = c(1,2), mar=c(1,1,1,1))

  VProfile(mod)

  if(withHist == TRUE) hist(rowSums(mod$V), col="darkred")

  VCProfile(mod)

  if(withHist == TRUE) hist(rowSums(mod$VC), col="darkred")
}


