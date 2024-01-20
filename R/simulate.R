

SIM = function(simObject, Tmax=200, TS=FALSE){
  states = saveStates(NULL, simObject)
  Lstates = saveLStates(NULL, simObject)

  for(i in 1:Tmax){
    simObject = AdultDynamics(simObject)
    simObject = AquaticDynamics(simObject)
    if(TS==TRUE) {
      states = saveStates(states, simObject)
      Lstates = saveLStates(Lstates, simObject)
    }
  }
  simObject$states = states
  simObject$Lstates = Lstates

  return(simObject)
}

steadyState = function(simObject, Tx=50, tol=.001){
  simObject = SIM(simObject, Tmax=Tx, TS=FALSE)

  err=10*tol
  while(err>tol){
    Bi = simObject$B; Qi = simObject$Q; Li = simObject$L
    simObject = SIM(simObject, Tmax=Tx, TS=FALSE)
    err = sum(abs(simObject$B - Bi)) +
      sum(abs(simObject$Q - Qi)) +
      sum(abs(simObject$L - Li))
  }

  steadyState = list()
  steadyState$B = simObject$B
  steadyState$Q = simObject$Q
  if(class(simObject) == "BQS") steadyState$S = simObject$S
  steadyState$L = simObject$L
  simObject$steadyState = steadyState
  return(simObject)
}




makeModel_BQ = function(b, q,
                        # Kernel Shapes, Search Weights
                        kFb, kFq,
                        wb=1, wq=1,
                        stayB=0, stayQ=0,
                        # Adult Parameters
                        adultPars = setup_BQ(),
                        # Aquatic Parameters
                        aquaticPars =  setupL(),
                        # Parasite Parameters
                        eip=15){

  model = makeSimObj_BQ(b, q, kFb, kFq, wb, wq, stayB, stayQ, adultPars, aquaticPars, eip)
  model = steadyState(model)
  model = makeKGV(model)
  model = makeTiles(model)
  model = makeAllGraphs(model)
  return(model)
}


makeModel_BQS = function(b, q, s,
                         # Kernel Shapes, Search Weights
                         kFb, kFq, kFs,
                         wb=1, wq=1, ws=1,
                         stayB=0, stayQ=0, stayS=0,
                         # Adult Parameters
                         adultPars = setup_BQS(),
                         # Aquatic Parameters
                         aquaticPars =  setupL(),
                         # Parasite Parameters
                         eip=15){

  model = makeSimObj_BQS(b, q, s, kFb, kFq, kFs, wb, wq, ws, stayB, stayQ, stayS, adultPars, aquaticPars, eip)
  model = steadyState(model)
  model = makeKGV(model)
  model = makeTiles(model)
  model = makeAllGraphs(model)
  return(model)
}

