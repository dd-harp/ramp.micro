
AdultDynamics.BQS = function(simObject){ with(simObject,{
  simObject$eggs = ova*psiQ*Q
  Bt = Mbb %*% B + Mqb %*% Q + Msb %*% S + Mlb %*% Lambda
  Qt = Mbq %*% B + Mqq %*% Q + Msq %*% S
  St = Mbs %*% B + Mqs %*% Q + Mss %*% S + Mls %*% Lambda
  browser()
  simObject$B = Bt
  simObject$Q = Qt
  simObject$S = St
  return(simObject)
})}

setup_BQS = function(pB=0.96, pQ=0.96, pS=0.96,
                     psiB=0.9, psiQ=0.9, psiS=0.9,
                     sigb=0.1, sigq=0.1, sigf=0.1, sigL=0.1, ova=20){
  par = list(pB=pB,pQ=pQ,pS=pS,
             psiB=psiB,psiQ=psiQ,psiS=psiS,
             sigb=sigb,sigq=sigq,sigf=sigf,sigL=sigL, ova=ova)
  class(par) <- "BQS"
  return(par)
}

init.BQS = function(simObject, B0=10, Q0=10, S0=10){
  simObject$B = matrix(B0, simObject$nb, 1)
  simObject$Q = matrix(Q0, simObject$nq, 1)
  simObject$S = matrix(S0, simObject$ns, 1)
  simObject$eggs = matrix(0, simObject$nq, 1)
  simObject = initL(simObject)
  return(simObject)
}

saveStates.BQS = function(states, simObject){
  if(is.null(states))
    return(list(Bt = simObject$B, Qt=simObject$Q, St=simObject$S))

  states = with(states,{
    Bt = cbind(Bt, simObject$B)
    Qt = cbind(Qt, simObject$Q)
    St = cbind(St, simObject$S)
    list(Bt=Bt, Qt=Qt, St=St)
  })

  return(states)
}


makeM_BQS = function(Pbb, Pqb, Psb,
                     Pbq, Pqq, Psq,
                     Pbs, Pqs, Pss,
                     pars){with(pars,{
                       nb = dim(Pbb)[1]
                       nq = dim(Pqq)[1]
                       ns = dim(Pss)[1]

                       Mbb = Pbb %*% diag(pB*(1-sigb)*(1-psiB), nb)
                       Mbq = Pbq %*% diag(pB*psiB, nb)
                       Mbs = Pbs %*% diag(pB*sigb*psiB, nb)

                       Mqb = Pqb %*% diag(pQ*(1-sigf)*psiQ, nq)
                       Mqq = Pqq %*% diag(pQ*(1-sigq)*(1-psiQ), nq)
                       Mqs = Pqs %*% diag(pQ*(sigf*psiQ + sigq*(1-psiQ)), nq)

                       Msb = Psb %*% diag(pS*psiS, ns)
                       Msq = 0*t(Mqs)
                       Mss = Pss %*% diag(pS*(1-psiS), ns)

                       Mlb = Pqb %*% diag(pQ*(1-sigL), nq)
                       Mls = Pqs %*% diag(pQ*sigL, nq)

                       bigM = rbind(
                         cbind(Mbb, Mqb, Msb),
                         cbind(Mbq, Mqq, Msq),
                         cbind(Mbs, Mqs, Mss)
                       )
                       return(list(Mbb=Mbb, Mbq=Mbq, Mbs=Mbs, Mqb=Mqb, Mqq=Mqq, Mqs=Mqs, Msb=Msb, Msq=Msq, Mss=Mss, Mlb=Mlb, Mls=Mls, bigM=bigM))
                     })}

makeSimObj_BQS = function(b, q, s, Kfb, Kfq, Kfs, wb=1, wq=1, ws=1,
                          stayB=0.1, stayQ=0.1, stayS=0.1,
                          # Adult Parameters
                          adultPars=setup_BQS(),
                          # Aquatic Parameters
                          aquaticPars = setupL(),
                          # Parasite Parameters
                          eip=15){

  Psi_bb = makePsi_stay(b,kFb,wb,stayB)
  Psi_bq = makePsi(b,q,kFq,wq)
  Psi_bs = makePsi(b,s,kFs,ws)
  Psi_qb = makePsi(q,b,kFb,wb)
  Psi_qq = makePsi_stay(q,kFq,wq,stayQ)
  Psi_qs = makePsi(q,s,kFs,ws)
  Psi_sb = makePsi(s,b,kFs,wb)
  Psi_sq = makePsi(s,q,kFs,wq)
  Psi_ss = makePsi_stay(s,kFs,ws,stayS)

  browser()
  simObject = makeM_BQS(Psi_bb, Psi_qb, Psi_sb, Psi_bq, Psi_qq, Psi_sq, Psi_bs, Psi_qs, Psi_ss, adultPars)

  simObject$b=b
  simObject$q=q
  simObject$s=s

  simObject$kFb = kFb
  simObject$kFq = kFq
  simObject$kFs = kFs

  nb = dim(b)[1]
  nq = dim(q)[1]
  ns = dim(s)[1]

  simObject$nb=nb
  simObject$nq=nq
  simObject$ns=ns

  simObject$Psi_bb = Psi_bb
  simObject$Psi_bq = Psi_bq
  simObject$Psi_bs = Psi_bs
  simObject$Psi_qb = Psi_qb
  simObject$Psi_qq = Psi_qq
  simObject$Psi_qs = Psi_qs
  simObject$Psi_sb = Psi_sb
  simObject$Psi_sq = Psi_sq
  simObject$Psi_ss = Psi_ss

  class(simObject) <- "BQS"

  simObject$adultPars= adultPars
  simObject$pB = adultPars$pB
  simObject$pQ = adultPars$pQ
  simObject$pS = adultPars$pS
  simObject$psiB = adultPars$psiB
  simObject$psiQ = adultPars$psiQ
  simObject$psiS = adultPars$psiS
  simObject$sigb = adultPars$sigb
  simObject$sigq = adultPars$sigq
  simObject$sigf = adultPars$sigf
  simObject$sigL = adultPars$sigL
  simObject$ova = adultPars$ova

  simObject$aquaticPars= aquaticPars
  simObject$pL = aquaticPars$pL
  simObject$theta = aquaticPars$theta
  simObject$zeta = aquaticPars$zeta
  simObject$xi = aquaticPars$xi

  simObject$eip = eip

  return(init(simObject))
}



makeKqb.BQS = function(simObject, Tmax=200){with(simObject,{
  nb = dim(Psi_bb)[2]
  nq = dim(Psi_qq)[2]
  ns = dim(Psi_ss)[2]

  Mbb = pB*(1-sigb)*(1-psiB)*Psi_bb
  success = pB*(1-sigb)*psiB*Psi_bb
  Mbq = pB*psiB*Psi_bq
  Mbs = pB*sigb*psiB*Psi_bs
  Mqb = pQ*(1-sigf)*psiQ*Psi_qb
  Mqq = pQ*(1-sigq)*(1-psiQ)*Psi_qq
  Mqs = pQ*(sigf*psiQ + sigq*(1-psiQ))*Psi_qs
  Msb = pS*psiS*Psi_sb
  Msq = 0*t(Mqs)
  Mss = (1-psiS)*Psi_ss

  M1 = cbind(Mbb, Mqb, Msb, 0*Mbb)
  M2 = cbind(0*Mbq, Mqq, Msq, 0*Mbq)
  M3 = cbind(Mbs, Mqs, Mss, 0*Mbs)
  M4 = cbind(success, 0*Mqb, 0*Msb, diag(1,nb))
  M = rbind(M1, M2, M3, M4)

  cohort = Psi_qb%*%diag(pQ,nq)
  Cyes = diag(psiB, nb)%*%cohort
  Cno = diag(1-psiB, nb)%*%cohort

  Kt = rbind(Cno, 0*Mqq, 0*Mqs, Cyes)

  for(i in 1:Tmax) Kt = M%*%Kt
  Kt[-c(1:(nb+nq+ns)),]

  simObject$Kqb = Kt[nb+nq+ns+c(1:nb),]
  return(simObject)
})}


makeKbq.BQS = function(simObject, Tmax = 200){with(simObject,{

  nb = dim(b)[1]
  nq = dim(q)[1]
  ns = dim(s)[1]

  Mbb = pB*(1-sigb)*(1-psiB)*Psi_bb
  success = pB*(1-sigb)*psiB*Psi_bb
  Mbq = pB*psiB*Psi_bq
  Mbs = pB*sigb*psiB*Psi_bs
  Mqb = pQ*(1-sigf)*psiQ*Psi_qb
  Mqq = pQ*(1-sigq)*(1-psiQ)*Psi_qq
  Mqs = pQ*(sigf*psiQ + sigq*(1-psiQ))*Psi_qs
  Msb = pS*psiS*Psi_sb
  Msq = 0*t(Mqs)
  Mss = (1-psiS)*Psi_ss

  M1 = cbind(Mbb, 0*Mqb, Msb, 0*Mqb)
  M2 = cbind(Mbq, Mqq, Msq, 0*Mqq)
  M3 = cbind(Mbs, Mqs, Mss, 0*Mqs)
  M4 = cbind(success, diag(1,nq), 0*Msq, diag(1,nq))
  M = rbind(M1, M2, M3, M4)

  cohort = Psi_bq%*%diag(pB,nb)
  Cyes = diag(psiQ, nq)%*%cohort
  Cno = diag(1-psiQ, nq)%*%cohort

  Kt = rbind(0*Mbb, Cno, 0*Mbs, Cyes)
  for(i in 1:Tmax) Kt = M%*%Kt
  Kt[-c(1:(nb+nq+ns)),]
  simObject$Kbq = Kt[nb+nq+ns+c(1:nq),]
  return(simObject)
})}

makeKbq_BQS = function(Pbb, Pqb, Psb,
                       Pbq, Pqq, Psq,
                       Pbs, Pqs, Pss,
                       sigf = 0.5, sigq = 0.5, sigb = 0.5,
                       pB =.98, pS = 0.99, pQ = 0.98,
                       psiB =.8, psiQ = 0.9, psiS = 0.98){
  nb = dim(b)[1]
  nq = dim(q)[1]
  ns = dim(s)[1]

  Mbb = pB*(1-sigb)*(1-psiB)*Pbb
  Mbq = pB*psiB*Pbq
  Mbs = pB*sigb*psiB*Pbs
  Mqb = pQ*(1-sigf)*psiQ*Pqb
  Mqq = pQ*(1-sigq)*(1-psiQ)*Pqq
  Mqs = pQ*(sigf*psiQ + sigq*(1-psiQ))*Pqs
  Msb = pS*psiS*Psb
  Msq = 0*t(Mqs)
  Mss = (1-psiS)*Pss

  M1 = cbind(Mbb, 0*Mqb, Msb, 0*Mqb)
  M2 = cbind(Mbq, Mqq, Msq, 0*Mqq)
  M3 = cbind(Mbs, 0*Mqs, Mss, 0*Mqs)
  M4 = cbind(0*Mbq, diag(1,nq), 0*Msq, diag(1,nq))
  M = rbind(M1, M2, M3, M4)

  Kt = rbind(0*Mbb, Mbq %*% diag(1,nb), 0*Mbs, 0*Mbq)
  for(i in 1:100) Kt = M%*%Kt
  Kt[-c(1:(nb+nq+ns)),]
}

