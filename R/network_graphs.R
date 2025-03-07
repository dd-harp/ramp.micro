# temporarily moved to RAMP-modelel-Library

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
  model$graphs$Kbb_net <- make_graph_obj(KGV$Kbb, "b", expression(K*scriptstyle(b%->%b)))
  model$graphs$Kqq_net <- make_graph_obj(KGV$Kqq, "q", expression(K*scriptstyle(q%->%q)))
  model$graphs$G_net   <- make_graph_obj(KGV$G,"q", "G")
  model$graphs$GG_net  <- make_graph_obj(KGV$GG,"q", "GG")
  model$graphs$V_net   <- make_graph_obj(KGV$V,"b", "V")
  model$graphs$VC_net  <- make_graph_obj(KGV$VC,"b", "VV")
  model$graphs$M_net   <- make_graph_obj(Mpar$bigM,"bq", "M")
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
  model$graphs$MM_net <-  make_graph_obj(bigMM,"bq", "MM")
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
  model$graphs$MM_net  <- make_graph_obj(bigMM, "bqs", "MM")
  return(model)
}


#' Plot the ouptuts of a graph using
#'
#' @param model a `ramp.micro` model object
#' @param net a network object
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
plot_graph = function(model, net, cut=NULL, alg="wt",
                      f_color = viridis::turbo,
                      min_edge_frac = 0.01, cx=2,
                      pw=1, mtl = "",
                      stretch=0.1,
                      r=.02, arw_lng=0.05, lwd=2){with(model,{
  if(alg == "wt") clusters = net$walktrap_clusters
  if(alg == "gr") clusters = net$greedy_clusters
  memix = if(is.null(cut)){
    membership(clusters)
  } else {
    cut_at(clusters, cut)
  }
  nC = max(memix)
  clrs = f_color(nC)
  frame_bq(b,q,mtl)
  if(net$type == "b"){
    xy = b
    add_arrows_xx(b, net$M, arw_clr = clrs[memix], min_edge_frac=min_edge_frac)
    add_points_b(b, wts_b=rowSums(net$M), clr=clrs[memix], cx_b=cx, pw=pw)
  }
  if(net$type == "q"){
    xy = q
    add_arrows_xx(q, net$M, arw_clr = clrs[memix], min_edge_frac=min_edge_frac)
    add_points_q(q, rowSums(net$M), clr=clrs[memix], cx_q=cx, pw=pw)
  }
  if(net$type == "q"){
    xy = rbind(b,q)
    add_arrows_xx(xy, net$M, arw_clr = clrs[memix], min_edge_frac=min_edge_frac)
    wts = rowSums(net$M)
    add_points_b(b, wts[1:nb], clr=clrs[memix[1:nb]], cx_b=cx, pw=pw)
    add_points_q(q, wts[nb+1:nq], clr=clrs[memix[nb+1:nq]], cx_q=cx, pw=pw)
  }
  add_convex_hulls(memix, xy, clrs, stretch, lwd)
  return(invisible())
})}
