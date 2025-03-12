
#' Plot point sets for a model
#'
#' @param model a compound [list] defining a model
#' @param wts_b values to scale the point size for points in b
#' @param wts_q values to scale the point size for points in q
#' @param wts_s values to scale the point size for points in s
#' @param cx_b set the maximum cex for b points
#' @param cx_q set the maximum cex for q points
#' @param cx_s set the maximum cex for s points
#' @param clr_b the color for b points
#' @param clr_q the color for q points
#' @param clr_s the color for s points
#' @param pw power relationship for scaling point size: pw=1 is linear
#' @param mtl a plot title
#'
#' @return invisible(NULL)
#' @export
plot_points = function(model,
                      wts_b=1, wts_q=1, wts_s=1,
                      cx_b=1, cx_q=1, cx_s=1,
                      clr_b="#fe5f55CC",
                      clr_q="skyblue",
                      clr_s="olivedrab2",
                      pw=1, mtl=""){
  UseMethod("plot_points", model$Mpar)
}

#' Plot point sets for a model for a BQ model
#'
#' @param model a compound [list] defining a model
#' @param wts_b values to scale the point size for points in b
#' @param wts_q values to scale the point size for points in q
#' @param wts_s values to scale the point size for points in s
#' @param cx_b set the maximum cex for b points
#' @param cx_q set the maximum cex for q points
#' @param cx_s set the maximum cex for s points
#' @param clr_b the color for b points
#' @param clr_q the color for q points
#' @param clr_s the color for s points
#' @param pw power relationship for scaling point size: pw=1 is linear
#' @param mtl a plot title
#'
#' @return invisible(NULL)
#' @export
plot_points.BQ = function(model,
                         wts_b=1, wts_q=1, wts_s=1,
                         cx_b=1, cx_q=1, cx_s=1,
                         clr_b="#fe5f55CC",
                         clr_q="skyblue",
                         clr_s="olivedrab2",
                         pw=1, mtl=""){
  with(model,plot_points_bq(b, q, wts_b, wts_q, cx_b, cx_q, clr_b, clr_q, pw, mtl))
  return(invisible())
}

#' Plot point sets for a bq model
#'
#' @param b a point set defining blood feeding sites
#' @param q a point set defining egg laying sites
#' @param wts_b values to scale the point size for points in b
#' @param wts_q values to scale the point size for points in q
#' @param cx_b set the maximum cex for b points
#' @param cx_q set the maximum cex for q points
#' @param clr_b the color for b points
#' @param clr_q the color for q points
#' @param pw power relationship for scaling point size: pw=1 is linear
#' @param mtl a plot title
#'
#' @return invisible(NULL)
#' @export
plot_points_bq = function(b, q,
                         wts_b=1, wts_q=1,
                         cx_b=1, cx_q=1,
                         clr_b="#fe5f55CC",
                         clr_q="skyblue",
                         pw=1,
                         mtl=""){
  frame_bq(b, q, mtl)
  add_points_b(b, wts_b, cx_b, clr_b, pw)
  add_points_q(q, wts_q, cx_q, clr_q, pw)
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
frame_bq = function(b, q, mtl=""){
  plot(rbind(b,q), type = "n", xaxt = "n", yaxt = "n",
       xlab = "", ylab = "", frame.plot=T, main = mtl)
  return(invisible())
}

#' Plot point sets for a model for a BQS model
#'
#' @param model a model defined as a compound [list]
#' @param wts_b values to scale the point size for points in b
#' @param wts_q values to scale the point size for points in q
#' @param wts_s values to scale the point size for points in s
#' @param cx_b set the maximum cex for b points
#' @param cx_q set the maximum cex for q points
#' @param cx_s set the maximum cex for s points
#' @param clr_b the color for b points
#' @param clr_q the color for q points
#' @param clr_s the color for s points
#' @param pw power relationship for scaling point size: pw=1 is linear
#' @param mtl a plot title
#'
#' @return no visible return value
#' @return invisible(NULL)
#' @export
plot_points.BQS = function(model,
                          wts_b=1, wts_q=1, wts_s=1,
                          cx_b=1, cx_q=1, cx_s=1,
                          clr_b="#fe5f55CC",
                          clr_q="skyblue",
                          clr_s="olivedrab2",
                          pw=1, mtl=""){
  with(model, plot_points_bqs(b, q, s,
                              wts_b, wts_q, wts_s,
                              cx_b, cx_q, cx_s,
                              clr_b, clr_q, clr_s,
                              pw, mtl))
  return(invisible())
}


#' Plot point sets for a model for a BQS model
#'
#' @param b a point set defining blood feeding sites
#' @param q a point set defining egg laying sites
#' @param s a point set defining sugar feeding sites
#' @param wts_b values to scale the point size for points in b
#' @param wts_q values to scale the point size for points in q
#' @param wts_s values to scale the point size for points in s
#' @param cx_b set the maximum cex for b points
#' @param cx_q set the maximum cex for q points
#' @param cx_s set the maximum cex for s points
#' @param clr_b the color for b points
#' @param clr_q the color for q points
#' @param clr_s the color for s points
#' @param pw power relationship for scaling point size: pw=1 is linear
#' @param mtl a plot title
#'
#' @return invisible(NULL)
#' @export
plot_points_bqs = function(b, q, s,
                           wts_b=1, wts_q=1, wts_s=1,
                           cx_b=1, cx_q=1, cx_s=1,
                           clr_b="#fe5f55CC",
                           clr_q="skyblue",
                           clr_s="olivedrab2",
                           pw=1, mtl=""){

  frame_bqs(b,q,s,mtl)
  add_points_s(s, wts_s, cx_s, clr_s, pw)
  add_points_b(b, wts_b, cx_b, clr_b, pw)
  add_points_q(q, wts_q, cx_q, clr_q, pw)
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
#' @param wts_b values to scale the point size
#' @param cx_b set the maximum cex for points
#' @param clr_b the color for b points
#' @param pw power relationship for scaling point size: pw=1 is linear
#'
#' @return invisible(NULL)
#' @export
add_points_b = function(b, wts_b=1, cx_b=1, clr_b="#fe5f55CC", pw=1){
  graphics::points(b, col = clr_b, pch = 15, cex=cx_b*wts_b^pw/max(wts_b^pw))
  return(invisible())
}

#' Add blood feeding points to a frame
#'
#' @param b a point set defining blood feeding sites
#' @param wts_b values to scale the point size
#' @param cx_b set the maximum cex for points
#' @param clr_b the outline color
#' @param pw power relationship for scaling point size: pw=1 is linear
#'
#' @return invisible(NULL)
#' @export
outline_points_b = function(b, wts_b=1, cx_b=1, clr_b="black", pw=1){
  graphics::points(b, pch=0, col=clr_b, cex=cx_b*wts_b^pw/max(wts_b^pw))
  return(invisible())
}

#' Add habitats to a frame
#'
#' @param q a point set defining egg laying sites
#' @param wts_q values to scale the point size
#' @param cx_q set the maximum cex for points
#' @param clr_q the color
#' @param pw power relationship for scaling point size: pw=1 is linear
#'
#' @return invisible(NULL)
#' @export
add_points_q = function(q, wts_q=1, cx_q=1, clr_q="skyblue", pw=1){
  graphics::points(q, col=clr_q, pch=19, cex=cx_q*wts_q^pw/max(wts_q^pw))
  return(invisible())
}

#' Add habitats to a frame
#'
#' @param q a point set defining egg laying sites
#' @param wts_q values to scale the point size
#' @param cx_q set the maximum cex for points
#' @param clr_q the color
#' @param pw power relationship for scaling point size: pw=1 is linear
#'
#' @return invisible(NULL)
#' @export
outline_points_q = function(q, wts_q=1, cx_q=1, clr_q="black", pw=1){
  graphics::points(q, col=clr_q, pch=1, cex=cx_q*wts_q^pw/max(wts_q^pw))
  return(invisible())
}

#' Add sugar feeding sites to a frame
#'
#' @param s a point set defining sugar feeding sites
#' @param wts_s values to scale the point size
#' @param cx_s set the maximum cex for points
#' @param clr_s the color
#' @param pw power relationship for scaling point size: pw=1 is linear
#'
#' @return invisible(NULL)
#' @export
add_points_s = function(s, wts_s=1, cx_s=1, clr_s="olivedrab2", pw=1){
  graphics::points(s, col=clr_s, pch=17, cex=cx_s*wts_s^pw/max(wts_s^pw))
  return(invisible())
}

#' Add habitats to a frame
#'
#' @param s a point set defining egg laying sites
#' @param wts_s values to scale the point size
#' @param cx_s set the maximum cex for points
#' @param clr_s the color
#' @param pw power relationship for scaling point size: pw=1 is linear
#'
#' @return invisible(NULL)
#' @export
outline_points_s = function(s, wts_s=1, cx_s=1, clr_s="black", pw=1){
  graphics::points(s, col=clr_s, pch=2, cex=cx_s*wts_s^pw/max(wts_s^pw))
  return(invisible())
}


#' Add blood feeding sites with and without the self loop
#'
#' @param b a point set defining blood feeding sites
#' @param M a matrix describing movement to b
#' @param cx_b set the maximum cex for points
#' @param clr_bA color with the self loop
#' @param clr_bB color without the self loop
#' @param pw power relationship for scaling point size: pw=1 is linear
#'
#' @return invisible(NULL)
#' @export
add_points_bb = function(b, M, cx_b=2, clr_bA="#cc444b66", clr_bB="#cc444bCC", pw=1){
  wts1 <- as.vector(rowSums(M))
  add_points_b(b, wts1, cx_b, clr_bA, pw)
  diag(M) <- 0
  wts2 <- as.vector(rowSums(M))
  add_points_b(b, wts2, cx_b*max(wts2)/max(wts1), clr_bB, pw)
  return(invisible())
}

#' Add habitats with and without the self loop
#'
#' @param q a point set defining egg laying sites
#' @param M a matrix describing movement to q
#' @param cx_q set the maximum cex for points
#' @param clr_qA color with the self loop
#' @param clr_qB color without the self loop
#' @param pw power relationship for scaling point size: pw=1 is linear
#'
#' @return invisible(NULL)
#' @export
add_points_qq = function(q, M, cx_q=2, clr_qA="skyblue", clr_qB="skyblue3", pw=1){
  wts1 <- as.vector(rowSums(M))
  add_points_q(q, wts1, cx_q, clr_qA, pw)
  diag(M) <- 0
  wts2 <- as.vector(rowSums(M))
  add_points_q(q, wts2, cx_q*max(wts2)/max(wts1), clr_qB, pw)
  return(invisible())
}

#' Add sugar feeding sites with and without the self loop
#'
#' @param s a point set defining sugar feeding sites
#' @param M a matrix describing movement to s
#' @param cx_s set the maximum cex for points
#' @param clr_sA color with the self loop
#' @param clr_sB color without the self loop
#' @param pw power relationship for scaling point size: pw=1 is linear
#'
#' @return invisible(NULL)
#' @export
add_points_ss = function(s, M, cx_s=2, clr_sA="olivedrab2", clr_sB="olivedrab3", pw=1){
  add_points_s(s, as.vector(rowSums(M)), cx_s, clr_sA, pw)
  diag(M) <- 0
  add_points_s(s, as.vector(rowSums(M)), cx_s, clr_sB, pw)
  return(invisible())
}




