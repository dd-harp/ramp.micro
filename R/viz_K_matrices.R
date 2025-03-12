
#' Plot the matrix \eqn{K_{q \leftarrow b}}
#' dispersal from \eqn{\left\{b\right\}} to \eqn{\left\{q\right\}}
#'
#' @param model a model defined as a compound [list]
#' @param cx_b set the maximum cex for b points
#' @param cx_q set the maximum cex for q points
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
plot_Kqb = function(model, cx_b=0.3, cx_q=2,
                    min_edge_frac=0.01, r=0.02, arw_lng=0.002, lwd=2,
                    clr_K="#4361eeCC", clr_b='red', clr_q ='darkblue'){
  with(model,with(Mpar,{
    frame_bq(b, q)
    add_arrows_xy(b, q, KGV$Kqb, min_edge_frac=min_edge_frac,
                  r=r, arw_lng=arw_lng, lwd=lwd, clr=clr_K)
    with(model, if(exists("s")) add_points_s(s, cx_s=cx_b))
    add_points_b(b, clr_b=clr_b, cx_b=cx_b)
    add_points_qq(q, KGV$Kqb, cx_q=cx_q, clr_qA=clr_q)
  }))
  return(invisible())
}

#' Plot the matrix \eqn{K_{b \leftarrow q}}: dispersal from \eqn{\left\{q\right\}} to \eqn{\left\{b\right\}}
#'
#' @param model a model defined as a compound [list]
#' @param cx_b set the maximum cex for b points
#' @param cx_q set the maximum cex for q points
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
plot_Kbq = function(model, cx_b=2, cx_q=0.3,
                    min_edge_frac=0.01, r=0.02, arw_lng=0.002, lwd=2,
                    clr_K="#fe5f55CC", clr_b='darkred', clr_q="#4361eeCC"){
  with(model,with(Mpar,{
    frame_bq(b, q)
    add_arrows_xy(q, b, KGV$Kbq, min_edge_frac=min_edge_frac,
                  r=r, arw_lng=arw_lng, lwd=lwd, clr=clr_K)
    with(model, if(exists("s")) add_points_s(s, cx_s=cx_q))
    add_points_q(q, cx_q=cx_q, clr_q=clr_q)
    add_points_bb(b, KGV$Kbq, cx_b=cx_b, clr_bB=clr_b)
  }))
  return(invisible())
}

#' Plot the matrix \eqn{K_{b \leftarrow b}}
#'
#' @param model a model defined as a compound [list]
#' @param cx_b set the maximum cex for b points
#' @param cx_q set the maximum cex for q points
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
plot_Kbb = function(model,  cx_b=2, cx_q=0.3,
                    min_edge_frac=0.01, r=0.02, arw_lng=0.002, lwd=2,
                    arw_clr ="#e2739655", seg_clr = '#00000022',
                    clr_b ="#cc444bCC", clr_q ="#4361eeCC"
){
  with(model,with(Mpar,
       {
         frame_bq(b, q)
         add_arrows_xx(b, KGV$Kbb, min_edge_frac=min_edge_frac,
                       r=r, arw_lng=arw_lng, lwd=lwd, arw_clr=arw_clr, seg_clr=seg_clr)
         with(model, if(exists("s")) add_points_s(s, cx_s=cx_q))
         add_points_q(q, cx_q=cx_q, clr_q=clr_q)
         add_points_bb(b, KGV$Kbb, cx_b=cx_b, clr_bA=arw_clr, clr_bB=clr_b)
       }))
  return(invisible())
}

#' Plot the matrix \eqn{K_{q \leftarrow q}}
#'
#' @param model a model defined as a compound [list]
#' @param cx_b set the maximum cex for b points
#' @param cx_q set the maximum cex for q points
#' @param min_edge_frac the minimum fraction to plot edge
#' @param r the radius of a ring around destination points
#' @param arw_lng the arrow length
#' @param lwd scale the line width
#' @param arw_clr the color to draw the arrow (asymetric part)
#' @param seg_clr the color to draw the segment (symmetric part)
#' @param clr_qA color for egg laying sites
#' @param clr_qB color for egg laying sites
#' @param clr_b color for blood feeding sites
#'
#' @return invisible(NULL)
#' @export
plot_Kqq = function(model, cx_b=0.03, cx_q=2,
                    min_edge_frac=0.01, r=0.02, arw_lng=0.002, lwd=2,
                    arw_clr = "#4361eeCC",
                    clr_qA= "#abc4ff55", seg_clr ='#00000022',
                    clr_qB="#858ae399", clr_b="#cc444bCC"){
  with(model,with(Mpar,{
    frame_bq(b, q)
    add_arrows_xx(q, KGV$Kqq, min_edge_frac=min_edge_frac,
                  r=r, arw_lng=arw_lng, lwd=lwd, arw_clr=arw_clr, seg_clr=seg_clr)
    with(model, if(exists("s")) add_points_s(s, cx_s = cx_b))
    add_points_b(b, cx_b = cx_b, clr_b=clr_b)
    add_points_qq(q, KGV$Kqq, cx_q=cx_q, clr_qA=clr_qA, clr_qB=clr_qB)
  }))
  return(invisible())
}
