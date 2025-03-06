# temporarily moved to RAMP-modelel-Library

#' Make graph object
#'
#' @description
#' From a matrix, M, return a graph object
#' the clusters from the walktrap algorithm &
#' for the symmetric graph defined by M + t(M)
#' the clusters for clusters_greedy
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
  clusters_walktrap <- cluster_walktrap(graph, weights=E(graph)$weight)
  MM = M+t(M)
  diag(MM) <-0
  symgraph <- graph_from_adjacency_matrix(MM, mode = "undirected", weighted=TRUE)
  clusters_greedy <- cluster_fast_greedy(symgraph, weights=E(symgraph)$weight)
  list(graph=graph, clusters_walktrap=clusters_walktrap, clusters_greedy=clusters_greedy, M=M, type = type, tag=tag)
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
