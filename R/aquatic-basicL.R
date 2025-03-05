
#' Aquatic Dynamics
#'
#' @inheritParams aquatic_dynamics
#'
#' @return the model, a compound [list]
#' @export
aquatic_dynamics.basicL = function(t, model){ with(model,
  with(c(Lpar, Lvars),{
    survive = pL*exp(-zeta*L)
    mature = theta*exp(-xi*L)
    L_t = (1-mature)*survive*L
    model$Lvars$L= L_t + model$terms$eggs
    model$terms$Lambda = mature*survive*L
    return(model)
}))}

#' Set initial values for the BQ model
#'
#' @inheritParams init_aquatic_model
#'
#' @return model, a compound [list]
#' @export
init_aquatic_model.basicL= function(model, L0_opts){
  model = init_aquatic_model_basicL(model, L0_opts)
  return(model)
}

#' Set initial values for
#'
#' @param model a model defined as a compound [list]
#' @param opts a list of options to overwrite defaults
#' @param L0 default initial value for larval density
#'
#' @return the model, a compound [list]
#' @export
init_aquatic_model_basicL = function(model, opts=list(), L0=10){with(opts,{
  model$Lvars$L = rep(L0, model$nq)[1:model$nq]
  model$terms$Lambda = rep(0, model$nq)[1:model$nq]
  return(model)
})}

#' Save state variables
#'
#' @inheritParams save_states_L
#'
#' @return the model, a compound [list]
#' @export
save_states_L.basicL = function(t, model){
  model$states$L$L_t[[t+1]] = model$Lvars$L
  model$states$L$Lambda_t[[t+1]] = model$Lvars$Lambda
  return(model)
}

#' The some of squared differences between two sets of variables
#'
#' @inheritParams compute_diffs_L
#'
#' @return a [numeric] value, the sum of squared differences
#' @export
compute_diffs_L.basicL = function(model){with(model,{
  dfs = sum((Lvars$L - steady$L$L)^2)
  return (dfs)
})}

#' Save the state variables in a vector
#'
#' @inheritParams init_states_L
#'
#' @return a [list] with the most recent states
#' @export
init_states_L.basicL = function(model){
  L = list()
  L$L_t = list()
  L$Lambda_t = list()
  model$states$L = L
  model = save_states_L(0, model)
  return(model)
}

#' Setup an aquatic model
#'
#' @inheritParams setup_aquatic_model
#'
#' @return a [list] defining an adult model
#' @export
setup_aquatic_model.basicL=function(model, aquatic_opts = list()){
  model = setup_aquatic_model_basicL(model, aquatic_opts)
  return(model)
}


#' Set up the aquatic population model of class basicL
#'
#' @param model a model defined as a compound [list]
#' @param opts a [list] of values that overwrite the defaults
#' @param pL base survival fraction
#' @param zeta reduced survival in response to mean crowding
#' @param theta base maturation fraction
#' @param xi delayed maturation in response to mean crowding
#'
#' @return the model, a compound [list]
#' @export
setup_aquatic_model_basicL = function(model, opts=list(),
                                      pL=.9, zeta=.01, theta=.9, xi=0){
  with(opts,{
    model$Lpar$pL=pL
    model$Lpar$zeta=zeta
    model$Lpar$theta=theta
    model$Lpar$xi=xi
    return(model)
})}
