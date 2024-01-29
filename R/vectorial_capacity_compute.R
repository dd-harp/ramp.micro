
#' Compute the potential dispersion of parasites after one feeding cycle
#'
#' @param model a compound [list] defining a model
#' @param Tmax the last time step
#'
#' @return the model, a compound [list]
#' @export
computeV = function(model, Tmax){
  UseMethod("computeV", model)
}

#' Compute the potential dispersion of parasites after one feeding cycle for the BQ model
#'
#' @param model a compound [list] defining a model
#' @param Tmax the last time step
#'
#' @return the model, a compound [list]
#' @export
computeV.BQ = function(model, Tmax=100){with(model,{
  Q = Mbq %*% diag(1, nb)
  B = diag(0, nb)

  for (i in 1:eip){
    Bt = Mbb %*% B + Mbq %*% Q
    Qt = Mqb %*% B + Mqq %*% Q
    B=Bt; Q=Qt
  }

  Vt = 0*B
  for (i in 1:Tmax){
    Vt = Vt + pB*psiB*B
    Bt = Mbb %*% B + Mbq %*% Q
    Qt = Mqb %*% B + Mqq %*% Q
    B=Bt; Q=Qt
  }
  model$V = Vt
  return(model)
})}

#' Compute the potential dispersion of parasites after one feeding cycle for the BQS model
#'
#' @param model a compound [list] defining a model
#' @param Tmax the last time step
#'
#' @return the model, a compound [list]
#' @export
computeV.BQS = function(model, Tmax=100){with(model,{
  Q = Mbq %*% diag(1, nb)
  B = diag(0, nb)
  S = matrix(0, ns, nb)

  for (i in 1:eip){
    Bt = Mbb %*% B + Mbq %*% Q + Mbs%*%S
    Qt = Mqb %*% B + Mqq %*% Q + Mqs%*%S
    St = Msb %*% B + Msq %*% Q + Mss%*%S
    B=Bt; Q=Qt
  }

  Vt = 0*B
  for (i in 1:Tmax){
    Vt = Vt + pB*psiB*B
    Bt = Mbb %*% B + Mbq %*% Q + Mbs%*%S
    Qt = Mqb %*% B + Mqq %*% Q + Mqs%*%S
    St = Msb %*% B + Msq %*% Q + Mss%*%S
    B=Bt; Q=Qt
  }
  model$V = Vt
  return(model)
})}

#' Compute the potential dispersion of parasites after one feeding cycle for the BQS model
#'
#' @param model a compound [list] defining a model
#'
#' @return the model, a compound [list]
#' @export
computeVC = function(model){with(model,{
  if(!exists("model$steadyState$B")) model = steadyState(model)
  model$VC = with(model, V %*% diag(as.vector(steadyState$B)))
  return(model)
})}


