
#' Compute the potential dispersion of parasites after one feeding cycle
#'
#' @param model a compound [list] defining a model
#' @param Tmax the last time step
#'
#' @return the model, a compound [list]
#' @export
compute_V = function(model, Tmax){
  UseMethod("compute_V", model$Mpar)
}

#' Compute the potential dispersion of parasites after one feeding cycle for the BQ model
#'
#' @param model a compound [list] defining a model
#' @param Tmax the last time step
#'
#' @return the model, a compound [list]
#' @export
compute_V.BQ = function(model, Tmax=100){with(model, with(Mpar,{
  Q = Mqb %*% diag(1, nb)
  B = diag(0, nb)

  for (i in 1:eip){
    Bt = Mbb %*% B + Mbq %*% Q
    Qt = Mqb %*% B + Mqq %*% Q
    B=Bt; Q=Qt
  }

  Vt = 0*B
  for (i in 1:Tmax){
    Vt = Vt + pB*psiB*B
    Bt = Mbb %*% B + Mbq %*% Q
    Qt = Mqb %*% B + Mqq %*% Q
    B=Bt; Q=Qt
  }
  model$KGV$V = Vt
  return(model)
}))}

#' Compute the potential dispersion of parasites after one feeding cycle for the BQS model
#'
#' @param model a compound [list] defining a model
#' @param Tmax the last time step
#'
#' @return the model, a compound [list]
#' @export
compute_V.BQS = function(model, Tmax=100){with(model, with(Mpar,{
  Q = Mqb %*% diag(1, nb)
  B = diag(0, nb)
  S = matrix(0, ns, nb)

  for (i in 1:eip){
    Bt = Mbb %*% B + Mbq %*% Q + Mbs%*%S
    Qt = Mqb %*% B + Mqq %*% Q + Mqs%*%S
    St = Msb %*% B + Msq %*% Q + Mss%*%S
    B=Bt; Q=Qt
  }

  Vt = 0*B
  for (i in 1:Tmax){
    Vt = Vt + pB*psiB*B
    Bt = Mbb %*% B + Mbq %*% Q + Mbs%*%S
    Qt = Mqb %*% B + Mqq %*% Q + Mqs%*%S
    St = Msb %*% B + Msq %*% Q + Mss%*%S
    B=Bt; Q=Qt
  }
  model$KGV$V = Vt
  return(model)
}))}

#' Compute the potential dispersion of parasites after one feeding cycle for the BQS model
#'
#' @param model a compound [list] defining a model
#'
#' @return the model, a compound [list]
#' @export
compute_VC = function(model){with(model,{
  if(!exists("model$steady$B")) model = steady_state(model)
  model$VC = with(model, Mpar$V %*% diag(as.vector(steady$M$B)))
  return(model)
})}


#' Visualize lifetime transmission from a site, per mosquito
#'
#' @param model a compound [list] defining a model
#' @param max_pt_sz_b the maximum cex for blood feeding sites
#' @param max_pt_sz_q the maximum cex for egg laying sites
#' @param min_edge_frac the fraction of the mass to plot
#' @param r the radius of a ring around destination points
#' @param arw_lng the arrow length
#' @param lwd scale the line width
#' @param lamp arrow width scaling factor
#' @param arw_clr the color to draw the arrow (asymetric part)
#' @param seg_clr the color to draw the segment (symmetric part)
#'
#' @return invisible(NULL)
#' @export
plot_dispersal_V = function(model,
                            max_pt_sz_b = 0.7, max_pt_sz_q = 2,
                            min_edge_frac = 0.01, r=.01, arw_lng=0.05, lwd=2, lamp=1,
                            arw_clr="darkolivegreen4", seg_clr="orangered3"){
  with(model,with(Mpar,{
    par(mar=c(2,2,2,2))
    frame_bq(b, q, mtl = "Potential Parasite Dispersal, per Mosquito")
    add_points_q(q, max_pt_sz = max_pt_sz_q)
    add_arrows_xx(b, KGV$V, min_edge_frac=min_edge_frac, r=r, arw_lng=arw_lng, lwd=lwd,
                  lamp=lamp, arw_clr=arw_clr, seg_clr=seg_clr)
    add_points_bb(b, KGV$V, max_pt_sz = max_pt_sz_b)
  }))
  return(invisible())
}

#' Visualize lifetime transmission from a site by mosquito populations
#'
#' @param model a compound [list] defining a model
#' @param max_pt_sz_b the maximum cex for blood feeding sites
#' @param max_pt_sz_q the maximum cex for egg laying sites
#' @param min_edge_frac the fraction of the mass to plot
#' @param r the radius of a ring around destination points
#' @param arw_lng the arrow length
#' @param lwd scale the line width
#' @param lamp arrow width scaling factor
#' @param arw_clr the color to draw the arrow (asymetric part)
#' @param seg_clr the color to draw the segment (symmetric part)
#'
#' @return invisible(NULL)
#' @export
plot_dispersal_VV = function(model,
                             max_pt_sz_b = 0.7, max_pt_sz_q = 2,
                             min_edge_frac = 0.01, r=.01, arw_lng=0.05, lwd=2, lamp=1,
                             arw_clr="springgreen4", seg_clr="firebrick3"){
  with(model,with(Mpar,{
    par(mar=c(2,2,2,2))
    frame_bq(b, q, mtl = "Potential Parasite Dispersal, Population")
    add_points_q(q, max_pt_sz = max_pt_sz_q)
    add_arrows_xx(b, KGV$VC, min_edge_frac=min_edge_frac, r=r, arw_lng=arw_lng, lwd=lwd,
                  lamp=lamp, arw_clr=arw_clr, seg_clr=seg_clr)
    add_points_bb(b, KGV$VC, max_pt_sz = max_pt_sz_b)
  }))
  return(invisible())
}


