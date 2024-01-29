
#' Plot the matrix Kqb: dispersal from {b} to {q}
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
  with(model,{
    frame_bq(b, q, mtl = expression(K[q%->%b]))
    if(exists("s")) add_points_s(s, max_pt_sz=0.5)
    add_arrows_xy(b, q, Kqb, min_edge_frac=min_edge_frac,
                  r=r, arw_lng=arw_lng, lwd=lwd, clr=clr_K)
    add_points_b(b, clr=clr_b, max_pt_sz=0.7)
    add_points_qq(q, Kqb, max_pt_sz=max_pt_sz, colA=clr_q)
  })
  return(invisible())
}

#' Plot the matrix Kbq: dispersal from {q} to {b}
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
  with(model,{
    frame_bq(b, q, mtl = expression(K[b%<-%q]))
    add_arrows_xy(b, q, Kbq, min_edge_frac=min_edge_frac,
                  r=r, arw_lng=arw_lng, lwd=lwd, clr=clr_K)
    with(model, if(exists("s")) add_points_s(s, max_pt_sz=0.7))
    add_points_q(q, max_pt_sz=0.7, clr=clr_q)
    add_points_bb(b, Kbq, pw=pw, max_pt_sz=max_pt_sz, colB=clr_b)
  })
  return(invisible())
}

#' Plot the matrix Kbb
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
  with(model,
       {
         frame_bq(b, q, mtl = expression(K[b %<-%b]))
         add_arrows_xx(b, Psi_bb, min_edge_frac=min_edge_frac,
                       r=r, arw_lng=arw_lng, lwd=lwd, arw_clr=arw_clr, seg_clr=seg_clr)
         with(model, if(exists("s")) add_points_s(s, max_pt_sz=0.7))
         add_points_q(q, max_pt_sz=0.7, clr=clr_q)
         add_points_bb(b, Kbb, max_pt_sz=max_pt_sz, colA=arw_clr, colB=clr_b)
       })
  return(invisible())
}

#' Plot the matrix Kqq
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
  with(model,{
    frame_bq(b, q, mtl = expression(K[q %<-%q]))
    with(model, if(exists("s")) add_points_s(s, max_pt_sz=0.7))
    add_arrows_xx(q, Psi_qq, min_edge_frac=min_edge_frac,
                  r=r, arw_lng=arw_lng, lwd=lwd, arw_clr=arw_clr, seg_clr=seg_clr)
    add_points_b(b, max_pt_sz=0.7, clr=clr_b)
    add_points_qq(q, Kqq, max_pt_sz=max_pt_sz2, colA=arw_clr, colB=clr_q)
  })
  return(invisible())
}
