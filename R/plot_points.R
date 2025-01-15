
#' Plot point sets for a model
#'
#' @param model a compound [list] defining a model
#' @param bwts values to scale the point size for points in b
#' @param qwts values to scale the point size for points in q
#' @param swts values to scale the point size for points in s
#' @param pw power relationship for scaling point size: pw=1 is linear
#' @param max_pt_sz set the maximum cex for points
#' @param mtl a plot title
#'
#' @return invisible(NULL)
#' @export
plot_points = function(model,
                      bwts=1, qwts=1, swts=1,
                      pw=1, max_pt_sz=1, mtl=""){
  UseMethod("plot_points", model$Mpar)
}

#' Plot point sets for a model for a BQ model
#'
#' @param model a compound [list] defining a model
#' @param bwts values to scale the point size for points in b
#' @param qwts values to scale the point size for points in q
#' @param swts values to scale the point size for points in s
#' @param pw power relationship for scaling point size: pw=1 is linear
#' @param max_pt_sz set the maximum cex for points
#' @param mtl a plot title
#'
#' @return invisible(NULL)
#' @export
plot_points.BQ = function(model,
                         bwts=1, qwts=1, swts=1,
                         pw=1, max_pt_sz=1, mtl=""){
  with(model,plot_points_bq(b, q, bwts, qwts, pw, max_pt_sz, mtl))
  return(invisible())
}

#' Plot point sets for a bq model
#'
#' @param b a point set defining blood feeding sites
#' @param q a point set defining egg laying sites
#' @param bwts values to scale the point size for points in b
#' @param qwts values to scale the point size for points in q
#' @param pw power relationship for scaling point size: pw=1 is linear
#' @param max_pt_sz set the maximum cex for points
#' @param mtl a plot title
#'
#' @return invisible(NULL)
#' @export
plot_points_bq = function(b, q,
                         bwts=1, qwts=1, pw=1, max_pt_sz=1,
                         mtl=""){
  frame_bq(b, q,mtl)
  add_points_b(b, bwts, pw, max_pt_sz)
  add_points_q(q, qwts, pw, max_pt_sz)
  return(invisible())
}

#' Set up a frame to plot points for a model
#'
#' @param b a point set defining blood feeding sites
#' @param q a point set defining egg laying sites
#' @param mtl a plot title
#'
#' @return invisible(NULL)
#' @export
frame_bq = function(b,q,mtl=""){
  plot(rbind(b,q), type = "n", xaxt = "n", yaxt = "n",
       xlab = "", ylab = "", frame.plot=T,
       main = mtl)
  return(NULL)
  return(invisible())
}

#' Plot point sets for a model for a BQS model
#'
#' @param model a model defined as a compound [list]
#' @param bwts values to scale the point size for points in b
#' @param qwts values to scale the point size for points in q
#' @param swts values to scale the point size for points in s
#' @param pw power relationship for scaling point size: pw=1 is linear
#' @param max_pt_sz set the maximum cex for points
#' @param mtl a plot title
#'
#' @return no visible return value
#' @return invisible(NULL)
#' @export
plot_points.BQS = function(model,
                          bwts=1, qwts=1, swts=1,
                          pw=1, max_pt_sz=1, mtl=""){
  with(model, plot_points_bqs(b, q, s, bwts, qwts, swts, pw, max_pt_sz, mtl))
  return(invisible())
}


#' Plot point sets for a model for a BQS model
#'
#' @param b a point set defining blood feeding sites
#' @param q a point set defining egg laying sites
#' @param s a point set defining sugar feeding sites
#' @param bwts values to scale the point size for points in b
#' @param qwts values to scale the point size for points in q
#' @param swts values to scale the point size for points in s
#' @param pw power relationship for scaling point size: pw=1 is linear
#' @param max_pt_sz set the maximum cex for points
#' @param mtl a plot title
#'
#' @return invisible(NULL)
#' @export
plot_points_bqs = function(b, q, s,
                          bwts=1, qwts=1, swts=1, pw=1, max_pt_sz=1,
                          mtl=""){

  frame_bqs(b,q,s,mtl)
  add_points_s(s, swts, pw, max_pt_sz)
  add_points_b(b, bwts, pw, max_pt_sz)
  add_points_q(q, qwts, pw, max_pt_sz)
  return(invisible())
}

#' Set up a frame to plot points for a BQS model
#'
#' @param b a point set defining blood feeding sites
#' @param q a point set defining egg laying sites
#' @param s a point set defining sugar feeding sites
#' @param mtl a plot title
#'
#' @return invisible(NULL)
#' @export
frame_bqs = function(b, q, s, mtl=""){
  plot(rbind(b,q,s), type = "n", xaxt = "n", yaxt = "n",
       xlab = "", ylab = "", frame.plot=T,
       main = mtl)
  return(invisible())
}

#' Add blood feeding points to a frame
#'
#' @param b a point set defining blood feeding sites
#' @param wts values to scale the point size
#' @param pw power relationship for scaling point size: pw=1 is linear
#' @param max_pt_sz set the maximum cex for points
#' @param clr the color
#'
#' @return invisible(NULL)
#' @export
add_points_b = function(b, wts=1, pw=1, max_pt_sz=1, clr="#fe5f55CC"){
  graphics::points(b, col = clr, pch = 15, cex=max_pt_sz*wts^pw/max(wts^pw))
  return(invisible())
}

#' Add habitats to a frame
#'
#' @param q a point set defining egg laying sites
#' @param wts values to scale the point size
#' @param pw power relationship for scaling point size: pw=1 is linear
#' @param max_pt_sz set the maximum cex for points
#' @param clr the color
#'
#' @return invisible(NULL)
#' @export
add_points_q = function(q, wts=1, pw=1, max_pt_sz=1, clr="skyblue"){
  graphics::points(q, col=clr, pch=19, cex=max_pt_sz*wts^pw/max(wts^pw))
  return(invisible())
}

#' Add sugar feeding sites to a frame
#'
#' @param s a point set defining sugar feeding sites
#' @param wts values to scale the point size
#' @param pw power relationship for scaling point size: pw=1 is linear
#' @param max_pt_sz set the maximum cex for points
#' @param clr the color
#'
#' @return invisible(NULL)
#' @export
add_points_s = function(s, wts=1, pw=1, max_pt_sz=1, clr="olivedrab2"){
  graphics::points(s, col=clr, pch=17, cex=max_pt_sz*wts^pw/max(wts^pw))
  return(invisible())
}

#' Add blood feeding sites with and without the self loop
#'
#' @param b a point set defining blood feeding sites
#' @param M a matrix describing movement to b
#' @param pw power relationship for scaling point size: pw=1 is linear
#' @param max_pt_sz set the maximum cex for points
#' @param colA color with the self loop
#' @param colB color without the self loop
#'
#' @return invisible(NULL)
#' @export
add_points_bb = function(b, M, pw=1, max_pt_sz=2, colA="#cc444b66", colB="#cc444bCC"){
  add_points_b(b, as.vector(rowSums(M)), pw, max_pt_sz, colA)
  diag(M) <- 0
  add_points_b(b, as.vector(rowSums(M)), pw, max_pt_sz, colB)
  return(invisible())
}

#' Add habitats with and without the self loop
#'
#' @param q a point set defining egg laying sites
#' @param M a matrix describing movement to q
#' @param pw power relationship for scaling point size: pw=1 is linear
#' @param max_pt_sz set the maximum cex for points
#' @param colA color with the self loop
#' @param colB color without the self loop
#'
#' @return invisible(NULL)
#' @export
add_points_qq = function(q, M, pw=1, max_pt_sz=2, colA="skyblue", colB="skyblue3"){
  add_points_q(q, as.vector(rowSums(M)), pw, max_pt_sz, colA)
  diag(M) <- 0
  add_points_q(q, as.vector(rowSums(M)), pw, max_pt_sz, colB)
  return(invisible())
}

#' Add sugar feeding sites with and without the self loop
#'
#' @param s a point set defining sugar feeding sites
#' @param M a matrix describing movement to s
#' @param pw power relationship for scaling point size: pw=1 is linear
#' @param max_pt_sz set the maximum cex for points
#' @param colA color with the self loop
#' @param colB color without the self loop
#'
#' @return invisible(NULL)
#' @export
add_points_ss = function(s, M, pw=1, max_pt_sz=2, colA="olivedrab2", colB="olivedrab3"){
  add_points_s(s, as.vector(rowSums(M)), pw, max_pt_sz, colA)
  diag(M) <- 0
  add_points_s(s, as.vector(rowSums(M)), pw, max_pt_sz, colB)
  return(invisible())
}




