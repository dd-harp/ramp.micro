
AquaticDynamics = function(simObject){
  UseMethod("AquaticDynamics", simObject$aquaticPars)
}

initL = function(simObject){
  UseMethod("initL", simObject$aquaticPars)
}

saveLStates=function(states, simObject){
  UseMethod("saveLStates", simObject$aquaticPars)
}

AquaticDynamics.L = function(simObject){ with(simObject,{
  survive = pL*exp(-zeta*L)
  mature = theta*exp(-xi*L)
  simObject$Lambda = mature*survive*L
  Lt = (1-mature)*survive*L
  Lt = Lt + eggs
  simObject$L= Lt
  return(simObject)
})}

setupL = function(theta=.9, zeta=.01, xi=0, pL=.9){
  p.L = list(theta=theta, zeta=zeta, xi=xi, pL=pL)
  class(p.L) <- "L"
  return(p.L)
}

initL.L = function(simObject, L0=10){
  simObject$L = matrix(L0, simObject$nq, 1)
  simObject$Lambda = matrix(0, simObject$nq, 1)
  return(simObject)
}

saveLStates.L = function(states, simObject){
  if(is.null(states))
    return(list(Lt = simObject$L))

  states = with(states,{
    Lt = cbind(Lt, simObject$L)
    list(Lt=Lt)
  })

  return(states)
}
