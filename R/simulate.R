
#' Simulate population dynamics
#'
#' @param model a compound [list] defining a model
#' @param Tmax the last time step
#'
#' @return model
#' @export
SIM = function(model, Tmax=200){
  model$states = list()
  model = init_states_M(model)
  model = init_states_L(model)

  for(t in 1:Tmax){
    model = AdultDynamics(t, model)
    model = AquaticDynamics(t, model)
    model = save_states_M(t, model)
    model = save_states_L(t, model)
  }

  return(model)
}

#' Run the model until it has reached a steady state
#'
#' @param model a compound [list] defining a model
#' @param burn initial burn time
#' @param Tx a run time chunk
#' @param tol tolerance for total sum of squared differences after a runtime chunk
#'
#' @return model
#' @export
steadyState = function(model, burn=500, Tx=50, tol=.001){

  for(t in 1:burn){
    model = AdultDynamics(t, model)
    model = AquaticDynamics(t, model)
  }

  steady=list()
  steady$M = model$Mvars
  steady$L = model$Lvars

  err=10*tol
  while(err>tol){
    for(t in 1:Tx){
      model = AdultDynamics(t, model)
      model = AquaticDynamics(t, model)
    }
    err = compute_diffs_M(model$Mvars, steady$M)
    err = err + compute_diffs_L(model$Lvars, steady$L)
    steady$M = model$Mvars
    steady$L = Model$Lvars
  }

  model$steady = steady

  return(model)
}
