
#' Compute dispersal matrix to lay eggs within one feeding cycle
#'
#' @param model a model defined as a compound [list]
#' @param Tmax the number of time steps
#'
#' @return the model, a compound [list]
#' @export
make_Kqb = function(model, Tmax){
  UseMethod("make_Kqb", model)
}

#' Compute dispersal matrix to blood feed within one feeding cycle
#'
#' @param model a model defined as a compound [list]
#' @param Tmax the number of time steps
#'
#' @return the model, a compound [list]
#' @export
make_Kbq = function(model, Tmax){
  UseMethod("make_Kbq", model)
}

#' Compute dispersal matrix to lay eggs within one feeding cycle
#'
#' @param model a model defined as a compound [list]
#'
#' @return the model, a compound [list]
#' @export
make_Kbb = function(model){
  model$Kbb = with(model, Kqb %*% Kbq)
  return(model)
}

#' Compute dispersal matrix to lay eggs within one feeding cycle
#'
#' @param model a model defined as a compound [list]
#'
#' @return the model, a compound [list]
#' @export
make_Kqq = function(model){
  model$Kqq = with(model, Kbq %*% Kqb)
  return(model)
}

#' Compute dispersal matrix to lay eggs within one feeding cycle
#'
#' @param model a model defined as a compound [list]
#' @param Tmax the number of time steps
#'
#' @return the model, a compound [list]
#' @export
make_Kqb.BQ = function(model, Tmax=100){with(model,{
  cohort = Psi_qb %*% diag(pQ, nq)
  Kqb = diag(psiB, nb) %*% cohort  # Success
  Bt = diag(1-psiB, nb) %*% cohort # Failure

  for(i in 1:Tmax){
    Kqb = Kqb + diag(pB*psiB, nb) %*% Psi_bb %*% Bt
    Bt = diag(pB*(1-psiB), nb) %*% Psi_bb %*% Bt
  }
  model$Kqb = Kqb
  return(model)
})}

#' Compute net dispersal matrix to blood feed within one feeding cycle
#'
#' @param model a model defined as a compound [list]
#' @param Tmax the number of time steps
#'
#' @return the model, a compound [list]
make_Kbq.BQ = function(model, Tmax=100){with(model,{
  cohort = Psi_bq %*% diag(pB, nb)
  Kbq = diag(psiQ, nq) %*% cohort  # Success
  Qt = diag(1-psiQ, nq) %*% cohort # Failure
  for(i in 1:Tmax){
    Kbq = Kbq + diag(pQ*psiQ, nq) %*% Psi_qq %*% Qt
    Qt =diag(pQ*(1-psiQ), nq) %*% Psi_qq %*%Qt
  }
  model$Kbq = Kbq
  return(model)
})}


#' Compute dispersal to lay eggs within one feeding cycle
#'
#' @param model a model defined as a compound [list]
#' @param Tmax the number of time steps
#'
#' @return the model, a compound [list]
#' @export
make_Kqb.BQS = function(model, Tmax=200){with(model,{
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

  model$Kqb = Kt[nb+nq+ns+c(1:nb),]
  return(model)
})}

#' Compute net dispersal matrix to blood feed within one feeding cycle
#'
#' @param model a model defined as a compound [list]
#' @param Tmax the number of time steps
#'
#' @return the model, a compound [list]
make_Kbq.BQS = function(model, Tmax = 200){with(model,{

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
  model$Kbq = Kt[nb+nq+ns+c(1:nq),]
  return(model)
})}


