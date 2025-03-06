
#' Make Psi from a source point set to a destination point set
#'
#' @param S a set of source points
#' @param D a set of destination points
#' @param kF a kernel weighting function
#' @param w linear weights for the destinations
#'
#' @return a [matrix], dispersal from S to D
#' @export
make_Psi_xy = function(S, D, kF=make_kF_exp(), w=1){
  lS = length(S[,1])
  lD = length(D[,1])
  K = matrix(0, lD, lS)
  if(length(w)==1) w=rep(w, lD)
  for(i in 1:lS){
    K[,i] = kF(sqrt((S[i,1]-D[,1])^2 + (S[i,2]-D[,2])^2), w)
    K[,i] = K[,i]/sum(K[,i])
  }
  return(K)
}

#' Make Psi from a source point set to a destination point set
#'
#' @param S a set of source points
#' @param kF a kernel weighting function
#' @param w linear weights for the destinations
#' @param stay the fraction that stays
#'
#' @return a [matrix], dispersal from S to S
#' @export
make_Psi_xx = function(S, kF=make_kF_exp(), w=1, stay=0){
  lS = dim(S)[1]
  K = matrix(0, lS, lS)
  if(length(w)==1) w=rep(w, lS)
  for(i in 1:lS){
    K[,i] = kF(sqrt((S[i,1]-S[,1])^2 + (S[i,2]-S[,2])^2), w)
    K[i,i] = 0
    K[,i] = (1-stay)*K[,i] /sum(K[,i])
    K[i,i] = stay
  }
  return(K)
}
