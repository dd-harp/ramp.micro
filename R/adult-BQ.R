
AdultDynamics.BQ = function(simObject){with(simObject,{
  simObject$eggs = ova*psiQ*Q
  Bt = Mbb %*% B + Mqb %*% Q + Mlb %*% Lambda
  Qt = Mbq %*% B + Mqq %*% Q
  simObject$B = Bt
  simObject$Q = Qt
  return(simObject)
})}

setup_BQ = function(pB=.96, pQ = 0.96, psiB=0.9, psiQ=0.9, ova=20){
  par = list(pB=pB,pQ=pQ,psiB=psiB,psiQ=psiQ,ova=ova)
  class(par) <- "BQ"
  return(par)
}

init.BQ = function(simObject, B0=10, Q0=10){
  simObject$B = matrix(B0, simObject$nb, 1)
  simObject$Q = matrix(Q0, simObject$nq, 1)
  simObject$eggs = matrix(0, simObject$nq, 1)
  simObject = initL(simObject)
  return(simObject)
}

saveStates.BQ = function(states, simObject){
  if(is.null(states))
    return(list(Bt = simObject$B, Qt = simObject$Q))

  states = with(states,{
    Bt = cbind(Bt, simObject$B)
    Qt = cbind(Qt, simObject$Q)
    list(Bt=Bt, Qt=Qt)
  })

  return(states)
}


makeM_BQ = function(Pbb, Pqb, Pbq, Pqq, pars){with(pars,{
  nb = dim(Pbb)[1]
  nq = dim(Pqq)[1]

  # "hardened" adults
  Mbb = Pbb %*% diag(pB*(1-psiB), nb)
  Mbq = Pbq %*% diag(pB*psiB, nb)
  Mqb = Pqb %*% diag(pQ*psiQ, nq)
  Mqq = Pqq %*% diag(pQ*(1-psiQ), nq)

  bigM = rbind(
    cbind(Mbb, Mqb),
    cbind(Mbq, Mqq)
  )
  # recently emerged adults
  Mlb = Pqb %*% diag(pQ, nq)

  list(Mbb=Mbb, Mbq=Mbq, Mqb=Mqb, Mqq=Mqq, Mlb=Mlb, bigM=bigM)
})}

makeSimObj_BQ = function(b, q,
                         # Kernel Shapes, Search Weights
                         kFb, kFq,
                         wb=1, wq=1,
                         stayB=0.1, stayQ=0.1,
                         # Adult Parameters
                         adultPars=setup_BQ(),
                         # Aquatic Parameters
                         aquaticPars=setupL(),
                         # Parasite Parameters
                         eip=15){
  Psi_bb = makePsi_stay(b,kFb,wb,stayB)
  Psi_bq = makePsi(b,q,kFq,wq)
  Psi_qb = makePsi(q,b,kFb,wb)
  Psi_qq = makePsi_stay(q,kFq,wq,stayQ)

  simObject = makeM_BQ(Psi_bb, Psi_qb, Psi_bq, Psi_qq, adultPars)

  simObject$b = b
  simObject$q = q
  simObject$kFb = kFb
  simObject$kFq = kFq
  simObject$nb = dim(b)[1]
  simObject$nq = dim(q)[1]

  simObject$Psi_bb = Psi_bb
  simObject$Psi_bq = Psi_bq
  simObject$Psi_qb = Psi_qb
  simObject$Psi_qq = Psi_qq

  class(simObject) <- "BQ"
  simObject$adultPars = adultPars
  simObject$pB = adultPars$pB
  simObject$pQ = adultPars$pQ
  simObject$psiB = adultPars$psiB
  simObject$psiQ = adultPars$psiQ
  simObject$ova = adultPars$ova

  simObject$aquaticPars = aquaticPars
  simObject$pL = aquaticPars$pL
  simObject$theta = aquaticPars$theta
  simObject$zeta = aquaticPars$zeta
  simObject$xi = aquaticPars$xi

  simObject$eip = eip

  return(init(simObject))
}



makeKqb.BQ = function(simObject, Tmax=100){with(simObject,{
  cohort = Psi_qb %*% diag(pQ, nq)
  Kqb = diag(psiB, nb) %*% cohort  # Success
  Bt = diag(1-psiB, nb) %*% cohort # Failure

  for(i in 1:Tmax){
    Kqb = Kqb + diag(pB*psiB, nb) %*% Psi_bb %*% Bt
    Bt = diag(pB*(1-psiB), nb) %*% Psi_bb %*% Bt
  }
  simObject$Kqb = Kqb
  return(simObject)
})}


makeKbq.BQ = function(simObject, Tmax=100){with(simObject,{
  cohort = Psi_bq %*% diag(pB, nb)
  Kbq = diag(psiQ, nq) %*% cohort  # Success
  Qt = diag(1-psiQ, nq) %*% cohort # Failure
  for(i in 1:Tmax){
    Kbq = Kbq + diag(pQ*psiQ, nq) %*% Psi_qq %*% Qt
    Qt =diag(pQ*(1-psiQ), nq) %*% Psi_qq %*%Qt
  }
  simObject$Kbq = Kbq
  return(simObject)
})}

