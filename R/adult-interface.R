

AdultDynamics = function(simObject){
  UseMethod("AdultDynamics", simObject)
}

init = function(simObject){
  UseMethod("init", simObject)
}

saveStates=function(states, simObject){
  UseMethod("saveStates", simObject)
}
makeKqb = function(simObject, Tmax){
  UseMethod("makeKqb", simObject)
}

makeKbq = function(simObject, Tmax){
  UseMethod("makeKbq", simObject)
}

