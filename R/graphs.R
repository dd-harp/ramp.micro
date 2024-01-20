

makeGraphs = function(M, type ="b", tag = ""){
  graph <- graph_from_adjacency_matrix(M, weighted=TRUE)
  clusters_walktrap <- cluster_walktrap(graph, weights=E(graph)$weight)
  MM = M+t(M)
  diag(MM) <-0
  symgraph <- graph_from_adjacency_matrix(MM, mode = "undirected", weighted=TRUE)
  clusters_greedy <- cluster_fast_greedy(symgraph, weights=E(symgraph)$weight)
  list(graph=graph, clusters_walktrap=clusters_walktrap,   clusters_greedy=clusters_greedy, M=M, type = type, tag=tag)
}

makeAllGraphs = function(mod){UseMethod("makeAllGraphs",mod)}

makeAllGraphs_common = function(mod){with(mod,{
  mod$Kbb_net <- makeGraphs(Kbb, "b", expression(K[b%->%b]))
  mod$Kqq_net <- makeGraphs(Kqq, "q", expression(K[q%->%q]))
  mod$G_net   <- makeGraphs(G,"q", "G")
  mod$GG_net  <- makeGraphs(GG,"q", "GG")
  mod$V_net   <- makeGraphs(V,"b", "V")
  mod$VC_net  <- makeGraphs(VC,"b", "VV")
  mod$M_net   <- makeGraphs(bigM,"bq", "M")
  return(mod)
})}

makeAllGraphs.BQ = function(mod){
  mod = makeAllGraphs_common(mod)
  mod$bigMM <- with(mod,bigM%*%diag(as.vector(c(steadyState$B, steadyState$Q))))
  mod$MM_net <-  with(mod,makeGraphs(bigMM,"bq", "MM"))
  return(mod)
}

makeAllGraphs.BQS = function(mod){
  mod = makeAllGraphs_common(mod)
  mod$bigMM <- with(mod,bigM%*%diag(as.vector(c(steadyState$B, steadyState$Q, steadyState$S))))
  mod$MM_net  <- with(mod,makeGraphs(bigMM, "bqs", "MM"))
  return(mod)
}
