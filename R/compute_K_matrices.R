
#' Compute dispersal matrix to lay eggs within one feeding cycle
#'
#' @param model a model defined as a compound [list]
#' @param Tmax the number of time steps
#'
#' @return the model, a compound [list]
#' @export
compute_Kqb = function(model, Tmax){
  UseMethod("compute_Kqb", model$Mpar)
}

#' Compute net dispersal matrix to blood feed within one feeding cycle
#'
#' @param model a model defined as a compound [list]
#' @param Tmax the number of time steps
#'
#' @return the model, a compound [list]
#' @export
compute_Kqb.BQ = function(model, Tmax=100){with(model, with(Mpar,{

  cohort = Psi_qb %*% diag(pB, nb)

  notpsiQ = diag(1-psiQ, nq)
  psiQ = diag(psiQ, nq)
  pQ = diag(pQ, nq)

  Kbq = 0*cohort  # Success
  Bt = cohort # Failure

  Kqb = psiQ %*% cohort  # Success
  Qt = notpsiQ %*% cohort # Failure
  for(i in 1:Tmax){
    Kqb = Kqb + psiQ %*% Qt
    Qt = Psi_qq %*% pQ %*% notpsiQ %*% Qt
  }
  model$KGV$Kqb = Kqb
  return(model)
}))}


#' Compute dispersal to lay eggs within one feeding cycle
#'
#' @param model a model defined as a compound [list]
#' @param Tmax the number of time steps
#'
#' @return the model, a compound [list]
#' @export
compute_Kqb.BQS = function(model, Tmax=200){with(model,with(Mpar,{

  nb = dim(b)[1]
  nq = dim(q)[1]
  ns = dim(s)[1]

  sigF = diag(sigf, nq)
  nsigF = diag(1-sigf, nq)
  npsiB = diag(1-psiB, nb)
  psiB = diag(psiB, nb)
  pB = diag(pB, nb)
  nsigB = diag(1-sigb, nb)
  sigB = diag(sigb, nb)

  npsiQ = diag(1-psiQ, nq)
  psiQ = diag(psiQ, nq)
  pQ = diag(pQ, nq)
  nsigQ = diag(1-sigq, nq)
  sigQ = diag(sigq, nq)

  npsiS = diag(1-psiS, ns)
  psiS = diag(psiS, ns)
  pS = diag(pS, ns)

  Mbb = Psi_bb %*% nsigB %*% pB %*% npsiB
  Mbq = 0*Psi_bq
  trap = psiQ
  Mbs = Psi_bs %*% psiS %*% pS

  Mqb = Psi_qb %*% psiB %*% pB
  Mqq = Psi_qq %*% nsigQ %*% pQ %*% npsiQ
  Mqs = 0*t(Psi_sq)

  Msb = Psi_sb %*% sigB %*% pB %*% npsiB
  Msq = Psi_sq %*% sigQ %*% pQ %*% npsiQ
  Mss = Psi_ss %*% npsiS %*% pS

  M1 = cbind(Mbb, Mbq, Mbs, 0*Mbq)
  M2 = cbind(Mqb, Mqq, Mqs, 0*Mqq)
  M3 = cbind(Msb, Msq, Mss, 0*Msq)
  M4 = cbind(0*Mqb, trap, 0*Mqs, diag(1,nq))
  M = rbind(M1, M2, M3, M4)

  B0 = pB
  Q0 = 0*(Psi_qb %*% pB)
  S0 = 0*(Psi_sb %*% pB)
  T0 = 0*Q0

  nix = nb + nq + ns + c(1:nq)
  BSQT = rbind(B0, S0, Q0, T0)

  looking = 1
  while(looking > 1e-6){
    BSQT = M%*%BSQT
    looking = sum(BSQT[-nix,-nix])
  }

  nix = 1:(nb+nq+ns)
  Kqb = BSQT[-nix,]
  model$KGV$Kqb = Kqb
  return(model)
}))}


#' Compute dispersal matrix to blood feed within one feeding cycle
#'
#' @param model a model defined as a compound [list]
#' @param Tmax the number of time steps
#'
#' @return the model, a compound [list]
#' @export
compute_Kbq = function(model, Tmax){
  UseMethod("compute_Kbq", model$Mpar)
}

#' Compute dispersal matrix to lay eggs within one feeding cycle
#'
#' @param model a model defined as a compound [list]
#' @param Tmax the number of time steps
#'
#' @return the model, a compound [list]
#' @export
compute_Kbq.BQ = function(model, Tmax=100){with(model, with(Mpar,{

  notpsiB = diag(1-psiB, nb)
  psiB = diag(psiB, nb)
  pB = diag(pB, nb)

  cohort = Psi_bq %*% diag(pQ, nq)
  Kbq = 0*cohort  # Success
  Bt = cohort # Failure

  for(i in 1:Tmax){
    Kbq = Kbq + psiB %*% Bt
    Bt = Psi_bb %*% pB %*% notpsiB %*% Bt
  }
  model$KGV$Kbq = Kbq
  return(model)
}))}


#' Compute net dispersal matrix to blood feed within one feeding cycle
#'
#' @param model a model defined as a compound [list]
#' @param Tmax the number of time steps
#'
#' @return the model, a compound [list]
#' @export
compute_Kbq.BQS = function(model, Tmax = 200){with(model,with(Mpar,{

  nb = dim(b)[1]
  nq = dim(q)[1]
  ns = dim(s)[1]

  sigF = diag(sigf, nq)
  nsigF = diag(1-sigf, nq)
  npsiB = diag(1-psiB, nb)
  psiB = diag(psiB, nb)
  pB = diag(pB, nb)
  nsigB = diag(1-sigb, nb)
  sigB = diag(sigb, nb)
  npsiS = diag(1-psiS, ns)
  psiS = diag(psiS, ns)
  pS = diag(pS, ns)

  cohort = diag(pQ,nq)
  B = Psi_bq %*% nsigF %*% diag(pQ,nq)
  S = Psi_sq %*% sigF %*% diag(pQ,nq)

  Mbb = Psi_bb %*% nsigB %*% pB %*% npsiB
  trap = psiB
  Msb = Psi_sb %*% sigB %*% pB %*% npsiB
  Mbs = Psi_bs %*% pS %*% psiS
  Mss = Psi_ss %*% pS %*% npsiS

  M1 = cbind(Mbb, Mbs, 0*Mbb)
  M2 = cbind(Msb, Mss, 0*Msb)
  M3 = cbind(trap, 0*Psi_bs, diag(1,nb))
  M = rbind(M1, M2, M3)

  nix = nb + ns + c(1:nb)
  TBS = rbind(B, S, 0*B)

  looking = 1
  while(looking > 1e-6) {
    TBS = M%*%TBS
    looking = sum(TBS[-nix,-nix])
  }
  nix = 1:(nb+ns)
  Kbq = TBS[-nix,]
  model$KGV$Kbq = Kbq
  return(model)
}))}


#' Compute dispersal matrix to lay eggs within one feeding cycle
#'
#' @param model a model defined as a compound [list]
#'
#' @return the model, a compound [list]
#' @export
compute_Kbb = function(model){
  model$KGV$Kbb = with(model$KGV, Kbq %*% Kqb)
  return(model)
}

#' Compute dispersal matrix to lay eggs within one feeding cycle
#'
#' @param model a model defined as a compound [list]
#'
#' @return the model, a compound [list]
#' @export
compute_Kqq = function(model){
  model$KGV$Kqq = with(model$KGV, Kqb %*% Kbq)
  return(model)
}



