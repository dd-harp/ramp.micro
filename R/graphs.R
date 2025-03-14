# temporarily moved to RAMP-modelel-Library

#' Choose one of the networks describing dispersion
#' @description
#' 1 -Kbb; 2 - Kqq;
#' 3 - G; 4 - GG;
#' 5 - V; 6 - VV;
#' 7 - M; 8 - MM
#'
#'
#' @param i an index for the type of model
#' @param model a model defined as a compound [list]
#'
#' @return a net object
#' @export
get_graph = function(model, i){
  if(i==1) graph = model$graphs$Kbb_graph
  if(i==2) graph = model$graphs$Kqq_graph
  if(i==3) graph = model$graphs$G_graph
  if(i==4) graph = model$graphs$GG_graph
  if(i==5) graph = model$graphs$V_graph
  if(i==6) graph = model$graphs$VC_graph
  if(i==7) graph = model$graphs$M_graph
  if(i==8) graph = model$graphs$MM_graph
  return(graph)
}


#' Make graph object
#'
#' @description
#' From a matrix, M, return a graph object
#' the clusters from the walktrap algorithm &
#' for the symmetric graph defined by M + t(M)
#' the clusters for greedy_clusters
#'
#'
#' @param M a matrix
#' @param type "b" for blood feeding; "q" for egg laying
#' @param tag a string for plotting
#'
#' @importFrom igraph graph_from_adjacency_matrix cluster_walktrap cluster_fast_greedy E
#' @returns a graph object
#' @export
make_graph_obj = function(M, type ="b", tag = ""){
  graph <- graph_from_adjacency_matrix(M, weighted=TRUE)
  walktrap_clusters <- cluster_walktrap(graph, weights=E(graph)$weight)
  MM = M+t(M)
  diag(MM) <-0
  symgraph <- graph_from_adjacency_matrix(MM, mode = "undirected", weighted=TRUE)
  greedy_clusters <- cluster_fast_greedy(symgraph, weights=E(symgraph)$weight)
  list(graph=graph, walktrap_clusters=walktrap_clusters, greedy_clusters=greedy_clusters, M=M, type = type, tag=tag)
}

#' Make all of the graphs
#'
#' @description
#'
#' Make common graphs for all the KGV objects and bigM
#'
#' @param model a ramp.micro model object
#'
#' @returns a ramp.micro model object
#' @export
make_common_graphs = function(model){with(model,{
  model$graphs <- list()
  model$graphs$Kbb_graph <- make_graph_obj(KGV$Kbb, "b", expression(K*scriptstyle(b%->%b)))
  model$graphs$Kqq_graph <- make_graph_obj(KGV$Kqq, "q", expression(K*scriptstyle(q%->%q)))
  model$graphs$G_graph   <- make_graph_obj(KGV$G,"q", "G")
  model$graphs$GG_graph  <- make_graph_obj(KGV$GG,"q", "GG")
  model$graphs$V_graph   <- make_graph_obj(KGV$V,"b", "V")
  model$graphs$VC_graph  <- make_graph_obj(KGV$VC,"b", "VV")
  model$graphs$M_graph   <- make_graph_obj(Mpar$bigM,"bq", "M")
  return(model)
})}

#' Make graphs
#'
#' @description
#' Make graphs for all the KGV objects and bigM
#'
#' @param model a ramp.micro model object
#'
#' @returns a ramp.micro model object
#' @export
make_all_graphs = function(model){
  UseMethod("make_all_graphs",model$Mpar)
}


#' Make graphs for a BQ model
#'
#' @description
#' Make graphs for all the KGV objects, bigM, and bigMM
#'
#' @param model a ramp.micro model object
#'
#' @returns a ramp.micro model object
#' @export
make_all_graphs.BQ = function(model){
  model = make_common_graphs(model)
  BQ <- with(model$steady$M, diag(as.vector(c(B, Q))))
  bigMM <- model$Mpar$bigM %*% BQ
  model$Mpar$bigMM <- bigMM
  model$graphs$MM_graphs <-  make_graph_obj(bigMM,"bq", "MM")
  return(model)
}

#' Make graphs for a BQS model
#'
#' @description
#' Make graphs for all the KGV objects and bigM
#'
#' @param model a ramp.micro model object
#'
#' @returns a ramp.micro model object
#' @export
make_all_graphs.BQS = function(model){
  model = make_common_graphs(model)
  BQS <- with(model$steady$M, diag(as.vector(c(B, Q, S))))
  bigMM <- model$Mpar$bigM %*% BQS
  model$Mpar$bigMM <- bigMM
  model$graphs$MM_graphs  <- make_graph_obj(bigMM, "bqs", "MM")
  return(model)
}


#' Plot the ouptuts of a graph using
#'
#' @param model a `ramp.micro` model object
#' @param graph a graph object
#' @param cut optional arguent for cut_at
#' @param alg walktrap = "wt" or greedy = "gr"
#' @param f_color a function that returns a list of colors (e.g. viridis::turbo)
#' @param min_edge_frac the fraction of the mass to plot
#' @param cx cex for plotting points
#' @param pw power relationship for scaling point size: pw=1 is linear
#' @param mtl a plot title
#' @param stretch make the hull slightly larger or slightly smaller
#' @param r the radius of a ring around destination points
#' @param arw_lng the arrow length
#' @param lwd scale the line width
#'
#' @return invisible(NULL)
#' @export
plot_graph = function(model, graph, cut=NULL, alg="wt",
                      f_color = viridis::turbo,
                      min_edge_frac = 0.01, cx=2,
                      pw=1, mtl = "",
                      stretch=0.1,
                      r=.02, arw_lng=0.05, lwd=2){with(model,{
  if(alg == "wt") clusters = graph$walktrap_clusters
  if(alg == "gr") clusters = graph$greedy_clusters
  memix = if(is.null(cut)){
    membership(clusters)
  } else {
    cut_at(clusters, cut)
  }
  nC = max(memix)
  clrs = f_color(nC)
  frame_bq(b,q,mtl)
  if(graph$type == "b"){
    xy = b
    add_arrows_xx(b, graph$M, arw_clr = clrs[memix], min_edge_frac=min_edge_frac)
    add_points_b(b, wts_b=rowSums(graph$M), clr=clrs[memix], cx_b=cx, pw=pw)
  }
  if(graph$type == "q"){
    xy = q
    add_arrows_xx(q, graph$M, arw_clr = clrs[memix], min_edge_frac=min_edge_frac)
    add_points_q(q, rowSums(graph$M), clr=clrs[memix], cx_q=cx, pw=pw)
  }
  if(graph$type == "bq"){
    xy = rbind(b,q)
    add_arrows_xx(xy, graph$M, arw_clr = clrs[memix], min_edge_frac=min_edge_frac)
    wts = rowSums(graph$M)
    add_points_b(b, wts[1:nb], clr=clrs[memix[1:nb]], cx_b=cx, pw=pw)
    add_points_q(q, wts[nb+1:nq], clr=clrs[memix[nb+1:nq]], cx_q=cx, pw=pw)
  }
#  add_convex_hulls(memix, xy, clrs, stretch, lwd)
  return(invisible())
})}

#' Plot the ouptuts of a graph using
#'
#' @param model a `ramp.micro` model object
#' @param graphs a graphswork object
#' @param cut optional arguent for cut_at
#' @param alg walktrap = "wt" or greedy = "gr"
#' @param f_color a function that returns a list of colors (e.g. viridis::turbo)
#' @param min_edge_frac the fraction of the mass to plot
#' @param cx cex for plotting points
#' @param pw power relationship for scaling point size: pw=1 is linear
#' @param mtl a plot title
#' @param stretch make the hull slightly larger or slightly smaller
#' @param r the radius of a ring around destination points
#' @param arw_lng the arrow length
#' @param lwd scale the line width
#'
#' @return invisible(NULL)
#' @export
plot_subgraph = function(model, graph, cut=NULL, alg="wt",
                      f_color = viridis::turbo,
                      min_edge_frac = 0.01, cx=2,
                      pw=1, mtl = "",
                      stretch=0.1,
                      r=.02, arw_lng=0.05, lwd=2){with(model,{
  if(alg == "wt") clusters = graph$walktrap_clusters
  if(alg == "gr") clusters = graph$greedy_clusters
  memix = if(is.null(cut)){
    membership(clusters)
  } else {
    cut_at(clusters, cut)
  }
  nC = max(memix)
  clrs = f_color(nC)
  frame_bq(b,q,mtl)
  for(i in 1:nC){
    mx <- max(rowSums(graph$M))
    ix = which(memix == i)
    M = graph$M[ix,][,ix]
    sx <- max(rowSums(M))
    if(graph$type == "b"){
      add_arrows_xx(b[ix,], M, arw_clr = clrs[i], lwd=lwd*sx/mx, min_edge_frac=min_edge_frac)
      add_points_b(b[ix,], wts_b=rowSums(M), clr=clrs[i], cx_b=cx*sx/mx, pw=pw)
      sxy = make_convex_hull_xy(b[ix,])
      sxy = stretch_convex_hull(sxy, 1+stretch)
      #graphics::polygon(sxy[,1], sxy[,2], border=clrs[i], lwd=lwd)
    }
    if(graph$type == "q"){
      add_arrows_xx(q[ix,], M, arw_clr = clrs[i], lwd=lwd*sx/mx, min_edge_frac=min_edge_frac)
      add_points_q(q[ix,], wts_q=rowSums(M), clr=clrs[i], cx_q=cx*sx/mx, pw=pw)
      sxy = make_convex_hull_xy(q[ix,])
      sxy = stretch_convex_hull(sxy, 1+stretch)
      #graphics::polygon(sxy[,1], sxy[,2], border=clrs[i], lwd=lwd)
    }
  }
  return(invisible())
                      })}

#' Plot the ouptuts of a graph using
#'
#' @param model a `ramp.micro` model object
#' @param graph a graph object
#' @param cut optional arguent for cut_at
#' @param alg walktrap = "wt" or greedy = "gr"
#' @param f_color a function that returns a list of colors (e.g. viridis::turbo)
#' @param min_edge_frac the fraction of the mass to plot
#' @param cx cex for plotting points
#' @param pw power relationship for scaling point size: pw=1 is linear
#' @param mtl a plot title
#' @param stretch make the hull slightly larger or slightly smaller
#' @param r the radius of a ring around destination points
#' @param arw_lng the arrow length
#' @param lwd scale the line width
#' @param lty the line tupe
#'
#' @return invisible(NULL)
#' @export
add_hulls = function(model, graph, cut=NULL, alg="wt",
                      f_color = viridis::turbo,
                      min_edge_frac = 0.01, cx=2,
                      pw=1, mtl = "",
                      stretch=0.1,
                      r=.02, arw_lng=0.05,
                      lwd=2, lty=1){with(model,{
  if(alg == "wt") clusters = graph$walktrap_clusters
  if(alg == "gr") clusters = graph$greedy_clusters
  memix = if(is.null(cut)){
    membership(clusters)
  } else {
    cut_at(clusters, cut)
  }
  nC = max(memix)
  clrs = f_color(nC)
  for(i in 1:nC){
    mx <- max(rowSums(graph$M))
    ix = which(memix == i)
    M = graph$M[ix,][,ix]
    sx <- max(rowSums(M))
    if(graph$type == "b"){
      sxy = make_convex_hull_xy(b[ix,])
      sxy = stretch_convex_hull(sxy, 1+stretch)
      graphics::polygon(sxy[,1], sxy[,2], border=clrs[i], lwd=lwd, lty=lty)
    }
    if(graph$type == "q"){
      sxy = make_convex_hull_xy(q[ix,])
      sxy = stretch_convex_hull(sxy, 1+stretch)
      graphics::polygon(sxy[,1], sxy[,2], border=clrs[i], lwd=lwd, lty=lty)
    }
  }
  return(invisible())
})}
