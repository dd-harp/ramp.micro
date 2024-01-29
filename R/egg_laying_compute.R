
#' Compute the dispersion of eggs after one feeding cycle
#'
#' @param model a compound [list] defining a model
#' @param Tmax the last time step
#'
#' @return the model, a compound [list]
#' @export
compute_G = function(model, Tmax){
  UseMethod("compute_G", model)
}

#' Compute the dispersion of eggs after one feeding cycle for the BQ model
#'
#' @param model a compound [list] defining a model
#' @param Tmax the last time step
#'
#' @return the model, a compound [list]
#' @export
compute_G.BQ = function(model, Tmax=50){with(model,{
  Q0 = diag(1, nq)
  B = Mqb %*% Q0
  Q = G = 0*Q0
  for(i in 1:Tmax){
    Bt = Mbb %*% B + Mqb %*% Q
    Qt = Mbq %*% B + Mqq %*% Q
    eggs = ova*psiQ*Q
    G = G + eggs
    B = Bt; Q=Qt
  }
  model$G = G
  return(model)
})}

#' Compute the dispersion of eggs after one feeding cycle for the BQ model
#'
#' @param model a compound [list] defining a model
#' @param Tmax the last time step
#'
#' @return the model, a compound [list]
#' @export
compute_G.BQS = function(model, Tmax=50){with(model,{
  Q0 = diag(1, nq)
  B = Mqb %*% Q0
  S = Mqs %*% Q0
  Q = G = 0*Q0
  for(i in 1:Tmax){
    Bt = Mbb %*% B + Mqb %*% Q + Msb%*%S
    Qt = Mbq %*% B + Mqq %*% Q + Msq%*%S
    St = Mbs %*% B + Mqs %*% Q + Mss%*%S
    eggs = ova*psiQ*Q
    G = G + eggs
    B = Bt; Q=Qt; S=St
  }
  model$G = G
  return(model)
})}

#' Compute the dispersion of eggs after one feeding cycle for the BQ model
#'
#' @param model a compound [list] defining a model
#'
#' @return the model, a compound [list]
#' @export
compute_GG = function(model){with(model,{
  if(!exists("model$steadyState$Q")) model = steadyState(model)
  model$GG = with(model,G %*% diag(as.vector(steadyState$Q)))
  return(model)
})}
