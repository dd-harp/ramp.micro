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
#' @param i an index for the type of model
#' @param model a model defined as a compound [list]
#'
#' @return a dispersion matrix
#' @export
getM_i = function(i, model){
  if(i==1) M = model$Kbb
  if(i==2) M = model$Kqq
  if(i==3) M = model$G
  if(i==4) M = model$GG
  if(i==5) M = model$V
  if(i==6) M = model$VC
  if(i==7) M = model$bigM
  if(i==8) M = model$bigMM
  return(M)
}

#' Choose one of the networks describing dispersion
#'
#' @param i an index for the type of model
#' @param model a model defined as a compound [list]
#'
#' @return a net object
#' @export
getNet.i = function(model, i){
  if(i==1) net = model$Kbb_net
  if(i==2) net = model$Kqq_net
  if(i==3) net = model$G_net
  if(i==4) net = model$GG_net
  if(i==5) net = model$V_net
  if(i==6) net = model$VC_net
  if(i==7) net = model$M_net
  if(i==8) net = model$MM_net
  return(net)
}
