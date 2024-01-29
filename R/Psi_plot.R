
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
plot_Psi = function(model, max_pt_sz=2,
                    min_edge_frac = 0.01, r=.01, arw_lng=0.05, lwd=2){
  UseMethod("plot_Psi", model)
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
plot_Psi.BQ = function(model, max_pt_sz=2,
                       min_edge_frac = 0.01, r=.01, arw_lng=0.05, lwd=2){
  with(model,{plot_Psi_BQ(b,q,Psi_bb, Psi_qb, Psi_bq, Psi_qq,
                          max_pt_sz=max_pt_sz, min_edge_frac=min_edge_frac,
                          r=r, arw_lng=arw_lng, lwd=lwd)})}

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
plot_Psi.BQS = function(model,max_pt_sz=2,
                        min_edge_frac = 0.01, r=.01, arw_lng=0.05, lwd=2){
  with(model,{plot_Psi_BQS(b,q,s,
               Psi_bb, Psi_qb, Psi_sb,
               Psi_bq, Psi_qq, Psi_sq,
               Psi_bs, Psi_qs, Psi_ss,
               max_pt_sz=max_pt_sz, min_edge_frac=min_edge_frac,
               r=r, arw_lng=arw_lng, lwd=lwd)
})}

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
plot_Psi_BQ = function(b,q,
                       Psi_bb, Psi_qb,
                       Psi_bq, Psi_qq,
                       max_pt_sz=2, min_edge_frac = 0.01,
                       r=.01, arw_lng=0.05, lwd=2){
  graphics::par(mfcol = c(2,2), mar = c(2,2,2,2))
  ## b to b
  frame_bq(b,q, mtl=expression(Psi[b %<-% b]))
  add_points_q(q, max_pt_sz = 0.7)
  add_arrows_xx(b, Psi_bb, min_edge_frac=min_edge_frac,
                r=r, arw_lng=arw_lng, lwd=lwd, arw_clr = "tomato")
  add_points_bb(b, Psi_bb, max_pt_sz=max_pt_sz)

  ## q to b
  frame_bq(b,q, mtl = expression(Psi[b %<-% q]))
  add_points_q(q, max_pt_sz = 0.7)
  add_arrows_xy(q, b, Psi_bq, min_edge_frac=min_edge_frac,
                r=r, arw_lng=arw_lng, lwd=lwd, clr="skyblue2")
  add_points_b(b, rowSums(Psi_bq), max_pt_sz = max_pt_sz)

  ## b to q
  frame_bq(b, q, mtl = expression(Psi[q %<-% b]))
  add_points_b(b, max_pt_sz=0.7)
  add_arrows_xy(b, q, Psi_qb, min_edge_frac=min_edge_frac,
                r=r, arw_lng=arw_lng, lwd=lwd, clr = "tomato")
  add_points_q(q, rowSums(Psi_qb), max_pt_sz=max_pt_sz)

  ##q to q
  frame_bq(b,q, mtl = expression(Psi[q %<-% q]))
  add_points_b(b, max_pt_sz=0.7)
  add_arrows_xx(q, Psi_qq, min_edge_frac=min_edge_frac,
                r=r, arw_lng=arw_lng, lwd=lwd, arw_clr="skyblue2")
  add_points_qq(q, Psi_qq, max_pt_sz=max_pt_sz, colA="skyblue")
  return(invisible())
}


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
plot_Psi_BQS = function(b,q,s,
                        Psi_bb, Psi_qb, Psi_sb,
                        Psi_bq, Psi_qq, Psi_sq,
                        Psi_bs, Psi_qs, Psi_ss,
                        max_pt_sz=2, min_edge_frac = 0.01,
                        r=.01, arw_lng=0.05, lwd=2){
  graphics::par(mfcol = c(3,3), mar = c(2,2,2,2))
  ## b to b
  frame_bqs(b,q,s, mtl = expression(Psi[b%<-%b]))
  add_arrows_xx(b, Psi_bb, min_edge_frac=min_edge_frac,
                r=r, arw_lng=arw_lng, lwd=lwd, arw_clr = "tomato")
  add_points_bb(b, Psi_bb, max_pt_sz=max_pt_sz)

  ## q to b
  frame_bqs(b, q, s, mtl = expression(Psi[b %<-% q]))
  add_points_q(q, max_pt_sz = 0.7)
  add_arrows_xy(q, b, Psi_bq, min_edge_frac=min_edge_frac,
                r=r, arw_lng=arw_lng, lwd=lwd, clr="skyblue2")
  add_points_b(b, rowSums(Psi_bq), max_pt_sz = max_pt_sz)

  ## s to b
  frame_bqs(b, q, s, mtl = expression(Psi[b %<-% s]))
  add_points_s(s, max_pt_sz = 0.7)
  add_arrows_xy(s, b, Psi_bs, min_edge_frac=min_edge_frac,
                r=r, arw_lng=arw_lng, lwd=lwd, clr="goldenrod")
  add_points_b(b, rowSums(Psi_bs), max_pt_sz = max_pt_sz)

  ## b to q
  frame_bqs(b, q, s, mtl = expression(Psi[q %<-% b]))
  add_points_b(b, max_pt_sz=0.7)
  add_arrows_xy(b, q, Psi_qb, min_edge_frac=min_edge_frac,
                r=r, arw_lng=arw_lng, lwd=lwd, clr = "tomato")
  add_points_q(q, rowSums(Psi_qb), max_pt_sz=max_pt_sz)

  ##q to q
  frame_bqs(b,q,s, mtl = expression(Psi[q%<-%q]))
  add_arrows_xx(q, Psi_qq, min_edge_frac=min_edge_frac,
                r=r, arw_lng=arw_lng, lwd=lwd, arw_clr = "skyblue2")
  add_points_qq(q, Psi_qq, max_pt_sz=max_pt_sz)

  ## s to q
  frame_bqs(b,q,s, mtl = expression(Psi[q %<-% s]))
  add_points_s(s, max_pt_sz = 0.7)
  add_arrows_xy(s, q, Psi_qs, min_edge_frac=min_edge_frac,
                r=r, arw_lng=arw_lng, lwd=lwd, clr="goldenrod")
  add_points_q(q, rowSums(Psi_qs), max_pt_sz = max_pt_sz)

  ## b to s
  frame_bqs(b, q, s, mtl = expression(Psi[s %<-% b]))
  add_points_b(b, max_pt_sz=0.7)
  add_arrows_xy(b, s, Psi_sb, min_edge_frac=min_edge_frac,
                r=r, arw_lng=arw_lng, lwd=lwd, clr = "tomato")
  add_points_s(s, rowSums(Psi_sb), max_pt_sz=max_pt_sz)

  ## q to s
  frame_bqs(b, q, s, mtl = expression(Psi[s %<-% q]))
  add_points_q(q, max_pt_sz=0.7)
  add_arrows_xy(q, s, Psi_sq, min_edge_frac=min_edge_frac,
                r=r, arw_lng=arw_lng, lwd=lwd, clr = "skyblue2")
  add_points_s(s, rowSums(Psi_sq), max_pt_sz=max_pt_sz)

  ##s to s
  frame_bqs(b,q,s, mtl = expression(Psi[s%<-%s]))
  add_arrows_xx(s, Psi_ss, min_edge_frac=min_edge_frac,
                r=r, arw_lng=arw_lng, lwd=lwd, arw_clr = "goldenrod")
  add_points_ss(s, Psi_ss, max_pt_sz=max_pt_sz)
  return(invisible())
}

