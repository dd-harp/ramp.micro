#' Decompose Matrix
#' @description Decompose a matrix into three parts: the
#' diagonal (self) loop, a symmetric matrix, and a flow
#'
#' @param M a matrix
#'
#' @return a [list] of the 3 decomposed matrices
#' @export
decompM = function(M){
  diagonal = diag(M)
  sym = pmin(M, t(M))
  flow = M-sym
  list(sym=sym, flow=flow, diagonal=diagonal)
}

#' Define a cutoff that includes a fraction of the total mass (e.g. 90 percent)
#'
#' @param M a matrix
#' @param fracMass the fraction of mass to plot
#'
#' @return a [numeric] value
#' @export
cutoffValue = function(M, fracMass=0.9){
  vals = sort(M, decreasing=TRUE)
  max(which(cumsum(vals)<fracMass*sum(M))) -> mx
  return(vals[mx])
}

#' Subset Edges
#' @description
#' Finds the indices of the edges that comprise a given
#' fraction of the total dispersal mass
#'
#' @param M a matrix
#' @param fracMass the fraction of mass to plot
#'
#' @return a set of array indices
#' @export
edgeSubset_fracMass = function(M, fracMass=0.9){
  V = cutoffValue(M, fracMass)
  return(which(M >= V, arr.ind=T))
}

#' Subset Edges
#' @description
#' Finds the indices of the edges that comprise a given
#' fraction of the total dispersal mass
#'
#' @param M a matrix
#' @param min_edge_frac the minimum fraction to plot edge
#'
#' @return a set of array indices
#' @export
edgeSubset = function(M, min_edge_frac=0.01){
  return(which(M >= min_edge_frac, arr.ind=T))
}

#' Choose one of the matrices describing movement
#'
#' @description
#' 1. Kbb; 2. Kqq;
#' 3. G; 4. GG;
#' 5. V; 6. VV;
#' 7. M; 8. MM;
#'
#' @param i an index for the type of model
#' @param model a model defined as a compound [list]
#'
#' @return a dispersion matrix
#' @export
get_matrix = function(i, model){
  if(i==1) M = model$KGV$Kbb
  if(i==2) M = model$KGV$Kqq
  if(i==3) M = model$KGV$G
  if(i==4) M = model$KGV$GG
  if(i==5) M = model$KGV$V
  if(i==6) M = model$KGV$VC
  if(i==7) M = model$Mpar$bigM
  if(i==8) M = model$Mpar$bigMM
  return(M)
}

#' Choose one of the networks describing dispersion
#' @description
#' 1. Kbb; 2. Kqq;
#' 3. G; 4. GG;
#' 5. V; 6. VV;
#' 7. M; 8. MM;
#'
#'
#' @param i an index for the type of model
#' @param model a model defined as a compound [list]
#'
#' @return a net object
#' @export
get_net = function(model, i){
  if(i==1) net = model$graphs$Kbb_net
  if(i==2) net = model$graphs$Kqq_net
  if(i==3) net = model$graphs$G_net
  if(i==4) net = model$graphs$GG_net
  if(i==5) net = model$graphs$V_net
  if(i==6) net = model$graphs$VC_net
  if(i==7) net = model$graphs$M_net
  if(i==8) net = model$graphs$MM_net
  return(net)
}
