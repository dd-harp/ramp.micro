
#' Aquatic Dynamics
#'
#' @param t the current time
#' @param model a model defined as a compound [list]
#'
#' @return the model, a compound [list]
#' @export
aquatic_dynamics = function(t, model){
  UseMethod("aquatic_dynamics", model$Lpar)
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
#' @param model a model
#'
#' @return a [numeric] value, the sum of squared differences
#' @export
compute_diffs_L = function(model){
  UseMethod("compute_diffs_L", model$Lpar)
}

#' Initialize the Lstates
#'
#' @param model a model defined as a compound [list]
#'
#' @return the model, a compound [list]
#' @export
init_states_L = function(model){
  UseMethod("init_states_L", model$Lpar)
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
