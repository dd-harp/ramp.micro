
#' Aquatic Dynamics
#'
#' @param t the current time
#' @param model a model defined as a compound [list]
#'
#' @return the model, a compound [list]
#' @export
AquaticDynamics = function(t, model){
  UseMethod("AquaticDynamics", model$Lpar)
}

#' Set initial values for
#'
#' @param model a model defined as a compound [list]
#' @param L0_opts a list of options to overwrite defaults
#'
#' @return the model, a compound [list]
#' @export
init_aquatic_model = function(model, L0_opts=list()){
  UseMethod("init_aquatic_model", model$Lpar)
}

#' Save state variables
#'
#' @param t the current time
#' @param model a model defined as a compound [list]
#'
#' @return the model, a compound [list]
#' @export
save_states_L=function(t, model){
  UseMethod("save_states_L", model$Lpar)
}

#' The some of squared differences between two sets of variables
#'
#' @param Lvars1 variables describing adult mosquitoes, set 1
#' @param Lvars2 variables describing adult mosquitoes, set 2
#'
#' @return a [numeric] value, the sum of squared differences
#' @export
compute_diffs_L = function(Lvars1, Lvars2){
  UseMethod("compute_diffs_L", Lvars1)
}

#' Initialize the Lstates
#'
#' @param model a model defined as a compound [list]
#'
#' @return the model, a compound [list]
#' @export
init_states_L = function(model){
  UseMethod("save_states_L", model$Lpar)
}

#' Setup an aquatic model
#'
#' @param model a model defined as a compound [list]
#' @param aquatic_opts a [list] to overwrite defaults
#'
#' @return a [list] defining an adult model
#' @export
setup_aquatic_model = function(model,aquatic_opts = list()){
  UseMethod("setup_aquatic_model", model$Lpar)
}
