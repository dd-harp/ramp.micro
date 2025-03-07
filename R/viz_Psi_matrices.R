
#' Visualize the one-bout dispersal matrices
#'
#' @param model a model defined as a compound [list]
#' @param cx_D the maximum cex for the destination
#' @param cx_S the maximum cex for the source
#' @param min_edge_frac the fraction of the mass to plot
#' @param r the radius of a ring around destination points
#' @param arw_lng the arrow length
#' @param lwd scale the line width
#'
#' @return the model, a compound [list]
#' @export
plot_all_Psi = function(model,
                        cx_D=2, cx_S=0.3,
                        min_edge_frac = 0.01,
                        r=.01, arw_lng=0.05, lwd=2){
  UseMethod("plot_all_Psi", model$Mpar)
}

#' Visualize the one-bout dispersal matrices for a BQ model
#'
#' @inheritParams plot_all_Psi
#'
#' @return the model, a compound [list]
#' @export
plot_all_Psi.BQ = function(model,
                           cx_D=2, cx_S=0.3,
                           min_edge_frac = 0.01,
                           r=.01, arw_lng=0.05, lwd=2){
  with(model,
    with(Mpar,{
      plot_all_Psi_BQ(b,q,
                      Psi_bb, Psi_qb, Psi_bq, Psi_qq,
                      cx_D=cx_D, cx_S=cx_S,
                      min_edge_frac=min_edge_frac,
                      r=r, arw_lng=arw_lng, lwd=lwd)
}))}



#' Visualize the one-bout dispersal matrices for a BQ model
#'
#' @param b blood feeding sites point set
#' @param q egg laying sites point set
#' @param Psi_bb one bout dispersal matrix from b to b
#' @param Psi_qb one bout dispersal matrix from q to b
#' @param Psi_bq one bout dispersal matrix from b to q
#' @param Psi_qq one bout dispersal matrix from q to q
#' @param cx_D the maximum cex for the destination
#' @param cx_S the maximum cex for the source
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
                       cx_D=2, cx_S=0.3,
                       min_edge_frac = 0.01,
                       r=.01, arw_lng=0.05, lwd=2){
  graphics::par(mfrow = c(2,2), mar = c(2,2,2,2))

  plot_Psi_bq(b, q, Psi_bq, cx_D, cx_S, min_edge_frac, r, arw_lng, lwd)
  plot_Psi_bb(b, q, Psi_bb, cx_D, min_edge_frac, r, arw_lng, lwd)
  plot_Psi_qb(b, q, Psi_qb, cx_D, cx_S, min_edge_frac, r, arw_lng, lwd)
  plot_Psi_qq(b, q, Psi_qq, cx_D, min_edge_frac, r, arw_lng, lwd)
  return(invisible())
}


#' Visualize the one-bout dispersal matrices for a BQS model
#'
#' @inheritParams plot_all_Psi
#'
#' @return the model, a compound [list]
#' @export
plot_all_Psi.BQS = function(model,
                            cx_D=2, cx_S=0.3,
                            min_edge_frac = 0.01,
                            r=.01, arw_lng=0.05, lwd=2){
  with(model,
     with(Mpar,{
       plot_all_Psi_BQS(b,q,s,
                        Psi_bb, Psi_qb, Psi_sb,
                        Psi_bq, Psi_qq, Psi_sq,
                        Psi_bs, Psi_qs, Psi_ss,
                        cx_D=cx_D, cx_S=cx_S,
                        min_edge_frac=min_edge_frac,
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
#' @param cx_D the maximum cex for the destination
#' @param cx_S the maximum cex for the source
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
                        cx_D=2, cx_S=0.3,
                        min_edge_frac = 0.01,
                        r=.01, arw_lng=0.05, lwd=2){
  graphics::par(mfcol = c(3,3), mar = c(2,2,2,2))

  plot_Psi_bb(b, q, Psi_bb, cx_D, min_edge_frac, r, arw_lng, lwd)
  plot_Psi_bq(b, q, Psi_bq, cx_D, cx_S, min_edge_frac, r, arw_lng, lwd)
  plot_Psi_bs(b, q, s, Psi_bs, cx_D, cx_S, min_edge_frac, r, arw_lng, lwd)
  plot_Psi_qb(b, q, Psi_qb, cx_D, cx_S, min_edge_frac, r, arw_lng, lwd)
  plot_Psi_qq(b, q, Psi_qq, cx_D, min_edge_frac, r, arw_lng, lwd)
  plot_Psi_qs(b, q, s, Psi_qs, cx_D, cx_S, min_edge_frac, r, arw_lng, lwd)
  plot_Psi_sb(b, q, s, Psi_sb, cx_D, cx_S, min_edge_frac, r, arw_lng, lwd)
  plot_Psi_sq(b, q, s, Psi_sq, cx_D, cx_S, min_edge_frac, r, arw_lng, lwd)
  plot_Psi_ss(b, q, s, Psi_ss, cx_D, min_edge_frac, r, arw_lng, lwd)
  return(invisible())
}

#' Visualize the one-bout dispersal matrices for a BQ model
#'
#' @param b blood feeding sites point set
#' @param q egg laying sites point set
#' @param Psi_bb one bout dispersal matrix from b to b
#' @param cx the maximum cex
#' @param min_edge_frac the fraction of the mass to plot
#' @param r the radius of a ring around destination points
#' @param arw_lng the arrow length
#' @param lwd scale the line width
#'
#' @return no visible return value
#' @export
plot_Psi_bb = function(b, q, Psi_bb,
                       cx=2, min_edge_frac = 0.01,
                       r=.02, arw_lng=0.05, lwd=2){

  ## b to b
  frame_bq(b, q, mtl=expression(Psi[b %<-% b]))
  add_arrows_xx(b, Psi_bb, min_edge_frac=min_edge_frac,
                r=r, arw_lng=arw_lng, lwd=lwd, arw_clr = "tomato")
  add_points_bb(b, Psi_bb, cx_b=cx)

  return(invisible())
}
#' Visualize the one-bout dispersal matrices for a BQ model
#'
#' @param b blood feeding sites point set
#' @param q egg laying sites point set
#' @param Psi_bq one bout dispersal matrix from b to q
#' @param cx_D the maximum cex for the destination
#' @param cx_S the maximum cex for the source
#' @param min_edge_frac the fraction of the mass to plot
#' @param r the radius of a ring around destination points
#' @param arw_lng the arrow length
#' @param lwd scale the line width
#'
#' @return no visible return value
#' @export
plot_Psi_bq = function(b, q, Psi_bq,
                       cx_D=2, cx_S=0.3,
                       min_edge_frac = 0.01,
                       r=.01, arw_lng=0.05, lwd=2){
  frame_bq(b, q, mtl = expression(Psi * scriptstyle(b %<-% q)))
  add_arrows_xy(q, b, Psi_bq, min_edge_frac=min_edge_frac,
                r=r, arw_lng=arw_lng, lwd=lwd, clr="tomato")
  add_points_q(q, cx_q = cx_S)
  add_points_b(b, rowSums(Psi_bq), cx_b = cx_D)
  return(invisible())
}

#' Visualize the dispersal matrix to b from s
#'
#' @param b blood feeding sites point set
#' @param q egg laying sites point set
#' @param s sugar feeding sites point set
#' @param Psi_bs one bout dispersal matrix to b from s
#' @param cx_D the maximum cex for the destination
#' @param cx_S the maximum cex for the source
#' @param min_edge_frac the fraction of the mass to plot
#' @param r the radius of a ring around destination points
#' @param arw_lng the arrow length
#' @param lwd scale the line width
#'
#' @return no visible return value
#' @export
plot_Psi_bs = function(b, q, s, Psi_bs,
                       cx_D=2, cx_S=0.3,
                       min_edge_frac = 0.01,
                       r=.01, arw_lng=0.05, lwd=2){
  frame_bqs(b, q, s, mtl = expression(Psi*scriptstyle(b %<-% s)))
  add_points_s(s, cx_s = cx_S)
  add_arrows_xy(s, b, Psi_bs, min_edge_frac=min_edge_frac,
                r=r, arw_lng=arw_lng, lwd=lwd, clr= "tomato")
  add_points_b(b, rowSums(Psi_bs), cx_b = cx_D)
  return(invisible())
}

#' Visualize the one-bout dispersal matrices for a BQ model
#'
#' @param b blood feeding sites point set
#' @param q egg laying sites point set
#' @param Psi_qb one bout dispersal matrix from q to b
#' @param cx_D the maximum cex for the destination
#' @param cx_S the maximum cex for the source
#' @param min_edge_frac the fraction of the mass to plot
#' @param r the radius of a ring around destination points
#' @param arw_lng the arrow length
#' @param lwd scale the line width
#'
#' @return no visible return value
#' @export
plot_Psi_qb = function(b, q, Psi_qb,
                       cx_D=2, cx_S=0.3,
                       min_edge_frac = 0.01,
                       r=.01, arw_lng=0.05, lwd=2){

  ## b to q
  frame_bq(b, q, mtl = expression(Psi*scriptstyle(q %<-% b)))
  add_arrows_xy(b, q, Psi_qb, min_edge_frac=min_edge_frac,
                r=r, arw_lng=arw_lng, lwd=lwd, clr = "skyblue")
  add_points_b(b, cx_b=cx_S)
  add_points_q(q, rowSums(Psi_qb), cx_q=cx_D)
  return(invisible())
}


#' Visualize the one-bout dispersal matrices for a BQ model
#'
#' @param b blood feeding sites point set
#' @param q egg laying sites point set
#' @param Psi_qq one bout dispersal matrix from q to q
#' @param cx the maximum cex
#' @param min_edge_frac the fraction of the mass to plot
#' @param r the radius of a ring around destination points
#' @param arw_lng the arrow length
#' @param lwd scale the line width
#'
#' @return no visible return value
#' @export
plot_Psi_qq = function(b, q, Psi_qq,
                       cx=2, min_edge_frac = 0.01,
                       r=.01, arw_lng=0.05, lwd=2){
  ##q to q
  frame_bq(b, q, mtl = expression(Psi*scriptstyle(q %<-% q)))
  add_arrows_xx(q, Psi_qq, min_edge_frac=min_edge_frac,
                r=r, arw_lng=arw_lng, lwd=lwd, arw_clr="skyblue")
  add_points_qq(q, Psi_qq, cx_q=cx)

  return(invisible())
}

#' Visualize the one-bout dispersal matrices for a BQ model
#'
#' @param b blood feeding sites point set
#' @param q egg laying sites point set
#' @param s sugar feeding sites point set
#' @param Psi_qs one bout dispersal matrix to q from s
#' @param cx_D the maximum cex for the destination
#' @param cx_S the maximum cex for the source
#' @param min_edge_frac the fraction of the mass to plot
#' @param r the radius of a ring around destination points
#' @param arw_lng the arrow length
#' @param lwd scale the line width
#'
#' @return no visible return value
#' @export
plot_Psi_qs = function(b, q, s, Psi_qs,
                       cx_D=2, cx_S=0.3,
                       min_edge_frac = 0.01,
                       r=.01, arw_lng=0.05, lwd=2){

  frame_bqs(b, q, s, mtl = expression(Psi*scriptstyle(q %<-% s)))
  add_arrows_xy(s, q, Psi_qs, min_edge_frac=min_edge_frac,
                r=r, arw_lng=arw_lng, lwd=lwd, clr = "skyblue")
  add_points_s(s, cx_s=cx_S)
  add_points_q(q, rowSums(Psi_qs), cx_q=cx_D)
  return(invisible())
}


#' Visualize the one-bout dispersal matrices for a BQ model
#'
#' @param b blood feeding sites point set
#' @param q egg laying sites point set
#' @param s sugar feeding sites point set
#' @param Psi_sb one bout dispersal matrix to s from b
#' @param cx_D the maximum cex for the destination
#' @param cx_S the maximum cex for the source
#' @param min_edge_frac the fraction of the mass to plot
#' @param r the radius of a ring around destination points
#' @param arw_lng the arrow length
#' @param lwd scale the line width
#'
#' @return no visible return value
#' @export
plot_Psi_sb = function(b, q, s, Psi_sb,
                       cx_D=2, cx_S=0.3,
                       min_edge_frac = 0.01,
                       r=.01, arw_lng=0.05, lwd=2){

  frame_bqs(b, q, s, mtl = expression(Psi*scriptstyle(s %<-% b)))
  add_arrows_xy(b, s, Psi_sb, min_edge_frac=min_edge_frac,
                r=r, arw_lng=arw_lng, lwd=lwd, clr = "olivedrab1")
  add_points_b(b, cx_b=cx_S)
  add_points_s(s, rowSums(Psi_sb), cx_s=cx_D)
  return(invisible())
}

#' Visualize the one-bout dispersal matrices for a BQ model
#'
#' @param b blood feeding sites point set
#' @param q egg laying sites point set
#' @param s sugar feeding sites point set
#' @param Psi_sq one bout dispersal matrix to s from q
#' @param cx_D the maximum cex for the destination
#' @param cx_S the maximum cex for the source
#' @param min_edge_frac the fraction of the mass to plot
#' @param r the radius of a ring around destination points
#' @param arw_lng the arrow length
#' @param lwd scale the line width
#'
#' @return no visible return value
#' @export
plot_Psi_sq = function(b, q, s, Psi_sq,
                       cx_D=2, cx_S=0.3,
                       min_edge_frac = 0.01,
                       r=.01, arw_lng=0.05, lwd=2){

  frame_bqs(b, q, s, mtl = expression(Psi*scriptstyle(s %<-% q)))
  add_arrows_xy(q, s, Psi_sq, min_edge_frac=min_edge_frac,
                r=r, arw_lng=arw_lng, lwd=lwd, clr = "olivedrab1")
  add_points_q(q, cx_q=cx_S)
  add_points_s(s, rowSums(Psi_sq), cx_s=cx_D)
  return(invisible())
}


#' Visualize the one-bout dispersal matrices for a Bs model
#'
#' @param b blood feeding sites point set
#' @param q egg laying sites point set
#' @param s sugar feeding sites point set
#' @param Psi_ss one bout dispersal matrix from s to s
#' @param cx the maximum cex
#' @param min_edge_frac the fraction of the mass to plot
#' @param r the radius of a ring around destination points
#' @param arw_lng the arrow length
#' @param lwd scale the line width
#'
#' @return no visible return value
#' @export
plot_Psi_ss = function(b, q, s, Psi_ss, cx=2,
                       min_edge_frac = 0.01,
                       r=.01, arw_lng=0.05, lwd=2){
  ##s to s
  frame_bqs(b, q, s, mtl = expression(Psi*scriptstyle(s %<-% s)))
  add_arrows_xx(s, Psi_ss, min_edge_frac=min_edge_frac,
                r=r, arw_lng=arw_lng, lwd=lwd, arw_clr="olivedrab1")
  add_points_ss(s, Psi_ss, cx_s=cx)

  return(invisible())
}
