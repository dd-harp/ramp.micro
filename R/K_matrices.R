
#' Compute dispersal matrix to lay eggs within one feeding cycle
#'
#' @param model a model defined as a compound [list]
#' @param Tmax the number of time steps
#'
#' @return the model, a compound [list]
#' @export
make_Kqb = function(model, Tmax){
  UseMethod("make_Kqb", model$Mpar)
}

#' Compute dispersal matrix to blood feed within one feeding cycle
#'
#' @param model a model defined as a compound [list]
#' @param Tmax the number of time steps
#'
#' @return the model, a compound [list]
#' @export
make_Kbq = function(model, Tmax){
  UseMethod("make_Kbq", model$Mpar)
}

#' Compute dispersal matrix to lay eggs within one feeding cycle
#'
#' @param model a model defined as a compound [list]
#'
#' @return the model, a compound [list]
#' @export
make_Kbb = function(model){
  model$KGV$Kbb = with(model$KGV, Kbq %*% Kqb)
  return(model)
}

#' Compute dispersal matrix to lay eggs within one feeding cycle
#'
#' @param model a model defined as a compound [list]
#'
#' @return the model, a compound [list]
#' @export
make_Kqq = function(model){
  model$KGV$Kqq = with(model$KGV, Kqb %*% Kbq)
  return(model)
}

#' Compute dispersal matrix to lay eggs within one feeding cycle
#'
#' @param model a model defined as a compound [list]
#' @param Tmax the number of time steps
#'
#' @return the model, a compound [list]
#' @export
make_Kbq.BQ = function(model, Tmax=100){with(model, with(Mpar,{
  cohort = Psi_bq %*% diag(pQ, nq)
  Kbq = diag(psiB, nb) %*% cohort  # Success
  Bt = diag(1-psiB, nb) %*% cohort # Failure

  for(i in 1:Tmax){
    Kbq = Kbq + diag(pB*psiB, nb) %*% Psi_bb %*% Bt
    Bt = diag(pB*(1-psiB), nb) %*% Psi_bb %*% Bt
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
make_Kqb.BQ = function(model, Tmax=100){with(model, with(Mpar,{
  cohort = Psi_qb %*% diag(pB, nb)
  Kqb = diag(psiQ, nq) %*% cohort  # Success
  Qt = diag(1-psiQ, nq) %*% cohort # Failure
  for(i in 1:Tmax){
    Kqb = Kqb + diag(pQ*psiQ, nq) %*% Psi_qq %*% Qt
    Qt =diag(pQ*(1-psiQ), nq) %*% Psi_qq %*%Qt
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
make_Kqb.BQS = function(model, Tmax=200){with(model,with(Mpar,{
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

  model$KGV$Kqb = Kt[nb+nq+ns+c(1:nb),]
  return(model)
}))}

#' Compute net dispersal matrix to blood feed within one feeding cycle
#'
#' @param model a model defined as a compound [list]
#' @param Tmax the number of time steps
#'
#' @return the model, a compound [list]
#' @export
make_Kbq.BQS = function(model, Tmax = 200){with(model,with(Mpar,{

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
  model$KGV$Kbq = Kt[nb+nq+ns+c(1:nq),]
  return(model)
}))}


#' Plot the matrix \eqn{K_{q \leftarrow b}}
#' dispersal from \eqn{\left\{b\right\}} to \eqn{\left\{q\right\}}
#'
#' @param model a model defined as a compound [list]
#' @param max_pt_sz set the maximum cex for points
#' @param min_edge_frac the fraction of the mass to plot
#' @param r the radius of a ring around destination points
#' @param arw_lng the arrow length
#' @param lwd scale the line width
#' @param clr_K the color for Kqb arrows
#' @param clr_b the color for b points
#' @param clr_q the color for q points
#'
#' @return invisible(NULL)
#' @export
plot_Kqb = function(model, max_pt_sz=2,
                    min_edge_frac=0.01, r=0.02, arw_lng=0.002, lwd=2,
                    clr_K="#4361eeCC", clr_b='red', clr_q ='darkblue'){
  with(model,with(Mpar,{
    frame_bq(b, q, mtl = expression(K*scriptstyle(q%<-%b)))
    add_arrows_xy(b, q, KGV$Kqb, min_edge_frac=min_edge_frac,
                  r=r, arw_lng=arw_lng, lwd=lwd, clr=clr_K)
    with(model, if(exists("s")) add_points_s(s, max_pt_sz=0.3))
    add_points_b(b, clr=clr_b, max_pt_sz=0.3)
    add_points_qq(q, KGV$Kqb, max_pt_sz=max_pt_sz, colA=clr_q)
  }))
  return(invisible())
}

#' Plot the matrix \eqn{K_{b \leftarrow q}}: dispersal from \eqn{\left\{q\right\}} to \eqn{\left\{b\right\}}
#'
#' @param model a model defined as a compound [list]
#' @param max_pt_sz set the maximum cex for points
#' @param min_edge_frac the fraction of the mass to plot
#' @param r the radius of a ring around destination points
#' @param arw_lng the arrow length
#' @param lwd scale the line width
#' @param clr_K the color for Kqb arrows
#' @param clr_b the color for b points
#' @param clr_q the color for q points
#'
#' @return invisible(NULL)
#' @export
plot_Kbq = function(model, max_pt_sz=2,
                    min_edge_frac=0.01, r=0.02, arw_lng=0.002, lwd=2,
                    clr_K="#fe5f55CC", clr_b='darkred', clr_q="#858ae399"){
  with(model,with(Mpar,{
    frame_bq(b, q, mtl = expression(K*scriptstyle(b%<-%q)))
    add_arrows_xy(q, b, KGV$Kbq, min_edge_frac=min_edge_frac,
                  r=r, arw_lng=arw_lng, lwd=lwd, clr=clr_K)
    with(model, if(exists("s")) add_points_s(s, max_pt_sz=0.3))
    add_points_q(q, max_pt_sz=0.3, clr=clr_q)
    add_points_bb(b, KGV$Kbq, max_pt_sz=max_pt_sz, colB=clr_b)
  }))
  return(invisible())
}

#' Plot the matrix \eqn{K_{b \leftarrow b}}
#'
#' @param model a model defined as a compound [list]
#' @param max_pt_sz set the maximum cex for points
#' @param min_edge_frac the fraction of the mass to plot
#' @param r the radius of a ring around destination points
#' @param arw_lng the arrow length
#' @param lwd scale the line width
#' @param arw_clr the color to draw the arrow (asymetric part)
#' @param seg_clr the color to draw the segment (symmetric part)
#' @param clr_b color for blood feeding sites
#' @param clr_q color for egg laying sites
#'
#' @export
plot_Kbb = function(model, max_pt_sz=2,
                    min_edge_frac=0.01, r=0.02, arw_lng=0.002, lwd=2,
                    arw_clr ="#e2739655", seg_clr = '#00000022',
                    clr_b ="#cc444bCC", clr_q ="#858ae399"
){
  with(model,with(Mpar,
       {
         frame_bq(b, q, mtl = expression(K*scriptstyle(b %<-%b)))
         add_arrows_xx(b, KGV$Kbb, min_edge_frac=min_edge_frac,
                       r=r, arw_lng=arw_lng, lwd=lwd, arw_clr=arw_clr, seg_clr=seg_clr)
         with(model, if(exists("s")) add_points_s(s, max_pt_sz=0.3))
         add_points_q(q, max_pt_sz=0.3, clr=clr_q)
         add_points_bb(b, KGV$Kbb, max_pt_sz=max_pt_sz, colA=arw_clr, colB=clr_b)
       }))
  return(invisible())
}

#' Plot the matrix \eqn{K_{q \leftarrow q}}
#'
#' @param model a model defined as a compound [list]
#' @param max_pt_sz set the maximum cex for points
#' @param min_edge_frac the minimum fraction to plot edge
#' @param r the radius of a ring around destination points
#' @param arw_lng the arrow length
#' @param lwd scale the line width
#' @param arw_clr the color to draw the arrow (asymetric part)
#' @param seg_clr the color to draw the segment (symmetric part)
#' @param clr_q color for egg laying sites
#' @param clr_b color for blood feeding sites
#'
#' @return invisible(NULL)
#' @export
plot_Kqq = function(model, max_pt_sz=2,
                    min_edge_frac=0.01, r=0.02, arw_lng=0.002, lwd=2,
                    arw_clr = "#abc4ff55", seg_clr ='#00000022',
                    clr_q="#858ae399", clr_b="#cc444bCC"){
  with(model,with(Mpar,{
    frame_bq(b, q, mtl = expression(K*scriptstyle(q %<-%q)))
    add_arrows_xx(q, KGV$Kqq, min_edge_frac=min_edge_frac,
                  r=r, arw_lng=arw_lng, lwd=lwd, arw_clr=arw_clr, seg_clr=seg_clr)
    with(model, if(exists("s")) add_points_s(s, max_pt_sz=0.7))
    add_points_b(b, max_pt_sz=0.3, clr=clr_b)
    add_points_qq(q, KGV$Kqq, max_pt_sz=max_pt_sz, colA=arw_clr, colB=clr_q)
  }))
  return(invisible())
}
