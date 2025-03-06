
#' Make Psi from a source point set to a destination point set
#'
#' @param S a set of source points
#' @param D a set of destination points
#' @param kF a kernel weighting function
#' @param w linear weights for the destinations
#'
#' @return a [matrix], dispersal from S to D
#' @export
make_Psi_xy = function(S, D, kF=make_kF_exp(), w=1){
  lS = length(S[,1])
  lD = length(D[,1])
  K = matrix(0, lD, lS)
  if(length(w)==1) w=rep(w, lD)
  for(i in 1:lS){
    K[,i] = kF(sqrt((S[i,1]-D[,1])^2 + (S[i,2]-D[,2])^2), w)
    K[,i] = K[,i]/sum(K[,i])
  }
  return(K)
}

#' Make Psi from a source point set to a destination point set
#'
#' @param S a set of source points
#' @param kF a kernel weighting function
#' @param w linear weights for the destinations
#' @param stay the fraction that stays
#'
#' @return a [matrix], dispersal from S to S
#' @export
make_Psi_xx = function(S, kF=make_kF_exp(), w=1, stay=0){
  lS = dim(S)[1]
  K = matrix(0, lS, lS)
  if(length(w)==1) w=rep(w, lS)
  for(i in 1:lS){
    K[,i] = kF(sqrt((S[i,1]-S[,1])^2 + (S[i,2]-S[,2])^2), w)
    K[i,i] = 0
    K[,i] = (1-stay)*K[,i] /sum(K[,i])
    K[i,i] = stay
  }
  return(K)
}

#' Visualize the one-bout dispersal matrices
#'
#' @param model a model defined as a compound [list]
#' @param max_pt_sz set the maximum cex for points
#' @param min_edge_frac the fraction of the mass to plot
#' @param r the radius of a ring around destination points
#' @param arw_lng the arrow length
#' @param lwd scale the line width
#'
#' @return the model, a compound [list]
#' @export
plot_all_Psi = function(model, max_pt_sz=2,
                    min_edge_frac = 0.01, r=.01, arw_lng=0.05, lwd=2){
  UseMethod("plot_all_Psi", model$Mpar)
}

#' Visualize the one-bout dispersal matrices for a BQ model
#'
#' @param model a model defined as a compound [list]
#' @param max_pt_sz set the maximum cex for points
#' @param min_edge_frac the fraction of the mass to plot
#' @param r the radius of a ring around destination points
#' @param arw_lng the arrow length
#' @param lwd scale the line width
#'
#' @return the model, a compound [list]
#' @export
plot_all_Psi.BQ = function(model, max_pt_sz=2,
                       min_edge_frac = 0.01, r=.01, arw_lng=0.05, lwd=2){
  with(model,with(Mpar,{plot_all_Psi_BQ(b,q,Psi_bb, Psi_qb, Psi_bq, Psi_qq,
                          max_pt_sz=max_pt_sz, min_edge_frac=min_edge_frac,
                          r=r, arw_lng=arw_lng, lwd=lwd)}))}



#' Visualize the one-bout dispersal matrices for a BQ model
#'
#' @param b blood feeding sites point set
#' @param q egg laying sites point set
#' @param Psi_bb one bout dispersal matrix from b to b
#' @param Psi_qb one bout dispersal matrix from q to b
#' @param Psi_bq one bout dispersal matrix from b to q
#' @param Psi_qq one bout dispersal matrix from q to q
#' @param max_pt_sz set the maximum cex for points
#' @param min_edge_frac the fraction of the mass to plot
#' @param r the radius of a ring around destination points
#' @param arw_lng the arrow length
#' @param lwd scale the line width
#'
#' @return no visible return value
#' @export
plot_all_Psi_BQ = function(b, q,
                       Psi_bb, Psi_qb,
                       Psi_bq, Psi_qq,
                       max_pt_sz=2, min_edge_frac = 0.01,
                       r=.01, arw_lng=0.05, lwd=2){
  graphics::par(mfrow = c(2,2), mar = c(2,2,2,2))

  plot_Psi_bq(b, q, Psi_bq, max_pt_sz, min_edge_frac, r, arw_lng, lwd)
  plot_Psi_bb(b, q, Psi_bb, max_pt_sz, min_edge_frac, r, arw_lng, lwd)
  plot_Psi_qb(b, q, Psi_qb, max_pt_sz, min_edge_frac, r, arw_lng, lwd)
  plot_Psi_qq(b, q, Psi_qq, max_pt_sz, min_edge_frac, r, arw_lng, lwd)
  return(invisible())
}


#' Visualize the one-bout dispersal matrices for a BQS model
#'
#' @param model a model defined as a compound [list]
#' @param max_pt_sz set the maximum cex for points
#' @param min_edge_frac the fraction of the mass to plot
#' @param r the radius of a ring around destination points
#' @param arw_lng the arrow length
#' @param lwd scale the line width
#'
#' @return the model, a compound [list]
#' @export
plot_all_Psi.BQS = function(model,max_pt_sz=2,
                        min_edge_frac = 0.01, r=.01, arw_lng=0.05, lwd=2){
  with(model,with(Mpar,{plot_all_Psi_BQS(b,q,s,
                           Psi_bb, Psi_qb, Psi_sb,
                           Psi_bq, Psi_qq, Psi_sq,
                           Psi_bs, Psi_qs, Psi_ss,
                           max_pt_sz=max_pt_sz, min_edge_frac=min_edge_frac,
                           r=r, arw_lng=arw_lng, lwd=lwd)
  }))}


#' Visualize the one-bout dispersal matrices for a BQS model
#'
#' @param b blood feeding sites point set
#' @param q egg laying sites point set
#' @param s sugar feeding sites point set
#' @param Psi_bb one bout dispersal matrix from b to b
#' @param Psi_qb one bout dispersal matrix from q to b
#' @param Psi_sb one bout dispersal matrix from s to b
#' @param Psi_bq one bout dispersal matrix from b to q
#' @param Psi_qq one bout dispersal matrix from q to q
#' @param Psi_sq one bout dispersal matrix from s to q
#' @param Psi_bs one bout dispersal matrix from b to s
#' @param Psi_qs one bout dispersal matrix from q to s
#' @param Psi_ss one bout dispersal matrix from s to s
#' @param max_pt_sz set the maximum cex for points
#' @param min_edge_frac the fraction of the mass to plot
#' @param r the radius of a ring around destination points
#' @param arw_lng the arrow length
#' @param lwd scale the line width
#'
#' @return no visible return value
#' @export
plot_all_Psi_BQS = function(b,q,s,
                        Psi_bb, Psi_qb, Psi_sb,
                        Psi_bq, Psi_qq, Psi_sq,
                        Psi_bs, Psi_qs, Psi_ss,
                        max_pt_sz=2, min_edge_frac = 0.01,
                        r=.01, arw_lng=0.05, lwd=2){
  graphics::par(mfcol = c(3,3), mar = c(2,2,2,2))

  plot_Psi_bb(b, q, Psi_bb, max_pt_sz, min_edge_frac, r, arw_lng, lwd)
  plot_Psi_bq(b, q, Psi_bq, max_pt_sz, min_edge_frac, r, arw_lng, lwd)
  plot_Psi_bs(b, q, s, Psi_bs, max_pt_sz, min_edge_frac, r, arw_lng, lwd)
  plot_Psi_qb(b, q, Psi_qb, max_pt_sz, min_edge_frac, r, arw_lng, lwd)
  plot_Psi_qq(b, q, Psi_qq, max_pt_sz, min_edge_frac, r, arw_lng, lwd)
  plot_Psi_qs(b, q, s, Psi_qs, max_pt_sz, min_edge_frac, r, arw_lng, lwd)
  plot_Psi_sb(b, q, s, Psi_sb, max_pt_sz, min_edge_frac, r, arw_lng, lwd)
  plot_Psi_sq(b, q, s, Psi_sq, max_pt_sz, min_edge_frac, r, arw_lng, lwd)
  plot_Psi_ss(b, q, s, Psi_ss, max_pt_sz, min_edge_frac, r, arw_lng, lwd)
  return(invisible())
}

#' Visualize the one-bout dispersal matrices for a BQ model
#'
#' @param b blood feeding sites point set
#' @param q egg laying sites point set
#' @param Psi_bb one bout dispersal matrix from b to b
#' @param max_pt_sz set the maximum cex for points
#' @param min_edge_frac the fraction of the mass to plot
#' @param r the radius of a ring around destination points
#' @param arw_lng the arrow length
#' @param lwd scale the line width
#'
#' @return no visible return value
#' @export
plot_Psi_bb = function(b, q, Psi_bb,
                       max_pt_sz=2, min_edge_frac = 0.01,
                       r=.02, arw_lng=0.05, lwd=2){

  ## b to b
  frame_bq(b,q, mtl=expression(Psi[b %<-% b]))
  add_arrows_xx(b, Psi_bb, min_edge_frac=min_edge_frac,
                r=r, arw_lng=arw_lng, lwd=lwd, arw_clr = "tomato")
  add_points_bb(b, Psi_bb, max_pt_sz=max_pt_sz)

  return(invisible())
}
#' Visualize the one-bout dispersal matrices for a BQ model
#'
#' @param b blood feeding sites point set
#' @param q egg laying sites point set
#' @param Psi_bq one bout dispersal matrix from b to q
#' @param max_pt_sz set the maximum cex for points
#' @param min_edge_frac the fraction of the mass to plot
#' @param r the radius of a ring around destination points
#' @param arw_lng the arrow length
#' @param lwd scale the line width
#'
#' @return no visible return value
#' @export
plot_Psi_bq = function(b,q, Psi_bq,
                       max_pt_sz=2, min_edge_frac = 0.01,
                       r=.01, arw_lng=0.05, lwd=2){
  frame_bq(b, q, mtl = expression(Psi * scriptstyle(b %<-% q)))
  add_arrows_xy(q, b, Psi_bq, min_edge_frac=min_edge_frac,
                r=r, arw_lng=arw_lng, lwd=lwd, clr="tomato")
  add_points_q(q, max_pt_sz = 0.3)
  add_points_b(b, rowSums(Psi_bq), max_pt_sz = max_pt_sz)
  return(invisible())
}

#' Visualize the dispersal matrix to b from s
#'
#' @param b blood feeding sites point set
#' @param q egg laying sites point set
#' @param s sugar feeding sites point set
#' @param Psi_bs one bout dispersal matrix to b from s
#' @param max_pt_sz set the maximum cex for points
#' @param min_edge_frac the fraction of the mass to plot
#' @param r the radius of a ring around destination points
#' @param arw_lng the arrow length
#' @param lwd scale the line width
#'
#' @return no visible return value
#' @export
plot_Psi_bs = function(b, q, s, Psi_bs,
                       max_pt_sz=2, min_edge_frac = 0.01,
                       r=.01, arw_lng=0.05, lwd=2){
  frame_bqs(b, q, s, mtl = expression(Psi*scriptstyle(b %<-% s)))
  add_points_s(s, max_pt_sz = 0.6)
  add_arrows_xy(s, b, Psi_bs, min_edge_frac=min_edge_frac,
                r=r, arw_lng=arw_lng, lwd=lwd, clr= "tomato")
  add_points_b(b, rowSums(Psi_bs), max_pt_sz = max_pt_sz)
  return(invisible())
}

#' Visualize the one-bout dispersal matrices for a BQ model
#'
#' @param b blood feeding sites point set
#' @param q egg laying sites point set
#' @param Psi_qb one bout dispersal matrix from q to b
#' @param max_pt_sz set the maximum cex for points
#' @param min_edge_frac the fraction of the mass to plot
#' @param r the radius of a ring around destination points
#' @param arw_lng the arrow length
#' @param lwd scale the line width
#'
#' @return no visible return value
#' @export
plot_Psi_qb = function(b, q, Psi_qb,
                       max_pt_sz=2, min_edge_frac = 0.01,
                       r=.01, arw_lng=0.05, lwd=2){

  ## b to q
  frame_bq(b, q, mtl = expression(Psi*scriptstyle(q %<-% b)))
  add_arrows_xy(b, q, Psi_qb, min_edge_frac=min_edge_frac,
                r=r, arw_lng=arw_lng, lwd=lwd, clr = "skyblue")
  add_points_b(b, max_pt_sz=0.3)
  add_points_q(q, rowSums(Psi_qb), max_pt_sz=max_pt_sz)
  return(invisible())
}


#' Visualize the one-bout dispersal matrices for a BQ model
#'
#' @param b blood feeding sites point set
#' @param q egg laying sites point set
#' @param Psi_qq one bout dispersal matrix from q to q
#' @param max_pt_sz set the maximum cex for points
#' @param min_edge_frac the fraction of the mass to plot
#' @param r the radius of a ring around destination points
#' @param arw_lng the arrow length
#' @param lwd scale the line width
#'
#' @return no visible return value
#' @export
plot_Psi_qq = function(b, q, Psi_qq,
                       max_pt_sz=2, min_edge_frac = 0.01,
                       r=.01, arw_lng=0.05, lwd=2){
  ##q to q
  frame_bq(b, q, mtl = expression(Psi*scriptstyle(q %<-% q)))
  add_arrows_xx(q, Psi_qq, min_edge_frac=min_edge_frac,
                r=r, arw_lng=arw_lng, lwd=lwd, arw_clr="skyblue")
  add_points_qq(q, Psi_qq, max_pt_sz=max_pt_sz)

  return(invisible())
}

#' Visualize the one-bout dispersal matrices for a BQ model
#'
#' @param b blood feeding sites point set
#' @param q egg laying sites point set
#' @param s sugar feeding sites point set
#' @param Psi_qs one bout dispersal matrix to q from s
#' @param max_pt_sz set the maximum cex for points
#' @param min_edge_frac the fraction of the mass to plot
#' @param r the radius of a ring around destination points
#' @param arw_lng the arrow length
#' @param lwd scale the line width
#'
#' @return no visible return value
#' @export
plot_Psi_qs = function(b, q, s, Psi_qs,
                       max_pt_sz=2, min_edge_frac = 0.01,
                       r=.01, arw_lng=0.05, lwd=2){

  frame_bqs(b, q, s, mtl = expression(Psi*scriptstyle(q %<-% s)))
  add_points_s(s, max_pt_sz=0.6)
  add_arrows_xy(s, q, Psi_qs, min_edge_frac=min_edge_frac,
                r=r, arw_lng=arw_lng, lwd=lwd, clr = "skyblue")
  add_points_q(q, rowSums(Psi_qs), max_pt_sz=max_pt_sz)
  return(invisible())
}


#' Visualize the one-bout dispersal matrices for a BQ model
#'
#' @param b blood feeding sites point set
#' @param q egg laying sites point set
#' @param s sugar feeding sites point set
#' @param Psi_sb one bout dispersal matrix to s from b
#' @param max_pt_sz set the maximum cex for points
#' @param min_edge_frac the fraction of the mass to plot
#' @param r the radius of a ring around destination points
#' @param arw_lng the arrow length
#' @param lwd scale the line width
#'
#' @return no visible return value
#' @export
plot_Psi_sb = function(b, q, s, Psi_sb,
                       max_pt_sz=2, min_edge_frac = 0.01,
                       r=.01, arw_lng=0.05, lwd=2){

  frame_bqs(b, q, s, mtl = expression(Psi*scriptstyle(s %<-% b)))
  add_points_b(b, max_pt_sz=0.6)
  add_arrows_xy(b, s, Psi_sb, min_edge_frac=min_edge_frac,
                r=r, arw_lng=arw_lng, lwd=lwd, clr = "olivedrab1")
  add_points_s(s, rowSums(Psi_sb), max_pt_sz=max_pt_sz)
  return(invisible())
}

#' Visualize the one-bout dispersal matrices for a BQ model
#'
#' @param b blood feeding sites point set
#' @param q egg laying sites point set
#' @param s sugar feeding sites point set
#' @param Psi_sq one bout dispersal matrix to s from q
#' @param max_pt_sz set the maximum cex for points
#' @param min_edge_frac the fraction of the mass to plot
#' @param r the radius of a ring around destination points
#' @param arw_lng the arrow length
#' @param lwd scale the line width
#'
#' @return no visible return value
#' @export
plot_Psi_sq = function(b, q, s, Psi_sq,
                       max_pt_sz=2, min_edge_frac = 0.01,
                       r=.01, arw_lng=0.05, lwd=2){

  frame_bqs(b, q, s, mtl = expression(Psi*scriptstyle(s %<-% q)))
  add_points_q(q, max_pt_sz=0.6)
  add_arrows_xy(q, s, Psi_sq, min_edge_frac=min_edge_frac,
                r=r, arw_lng=arw_lng, lwd=lwd, clr = "olivedrab1")
  add_points_s(s, rowSums(Psi_sq), max_pt_sz=max_pt_sz)
  return(invisible())
}


#' Visualize the one-bout dispersal matrices for a Bs model
#'
#' @param b blood feeding sites point set
#' @param q egg laying sites point set
#' @param s sugar feeding sites point set
#' @param Psi_ss one bout dispersal matrix from s to s
#' @param max_pt_sz set the maximum cex for points
#' @param min_edge_frac the fraction of the mass to plot
#' @param r the radius of a ring around destination points
#' @param arw_lng the arrow length
#' @param lwd scale the line width
#'
#' @return no visible return value
#' @export
plot_Psi_ss = function(b, q, s, Psi_ss,
                       max_pt_sz=2, min_edge_frac = 0.01,
                       r=.01, arw_lng=0.05, lwd=2){
  ##s to s
  frame_bqs(b, q, s, mtl = expression(Psi*scriptstyle(s %<-% s)))
  add_arrows_xx(s, Psi_ss, min_edge_frac=min_edge_frac,
                r=r, arw_lng=arw_lng, lwd=lwd, arw_clr="olivedrab1")
  add_points_ss(s, Psi_ss, max_pt_sz=max_pt_sz)

  return(invisible())
}

