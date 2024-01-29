#' Adult dynamics
#'
#' @param model a model defined as a compound [list]
#'
#' @param t the current time
#' @return the model, a compound [list]
#' @export
adult_dynamics = function(t, model){
  UseMethod("adult_dynamics", model$Mpar)
}

#' Save state variables
#'
#' @param t the current time
#' @param model a model defined as a compound [list]
#'
#' @return the model, a compound [list]
#' @export
save_states_M = function(t, model){
  UseMethod("save_states_M", model$Mpar)
}

#' Save state variables
#'
#' @param model a model defined as a compound [list]
#'
#' @return the model, a compound [list]
#' @export
init_states_M = function(model){
  UseMethod("init_states_M", model$Mpar)
}


#' Set initial values for
#'
#' @param model a model defined as a compound [list]
#' @param M0_opts a list of options to overwrite defaults
#'
#' @return the model, a compound [list]
#' @export
init_adult_model = function(model, M0_opts){
  UseMethod("init_adult_model", model$Mpar)
}

#' The some of squared differences between two sets of variables
#'
#' @param Mvars1 variables describing adult mosquitoes, set 1
#' @param Mvars2 variables describing adult mosquitoes, set 2
#'
#' @return a [numeric] value, the sum of squared differences
#' @export
compute_diffs_M = function(Mvars1, Mvars2){
  UseMethod("compute_diffs_M", Mvars1)
}

#' Setup an adult model
#'
#' @param model a model defined as a compound [list]
#' @param b a point set defining blood feeding sites
#' @param q a point set defining egg laying sites
#' @param s a point set defining sugar feeding sites
#' @param dispersal_opts a [list] to overwrite defaults
#' @param bionomic_opts a [list] to overwrite defaults
#' @param eip the extrinsic incubation period
#'
#' @return a [list] defining an adult model
#' @export
setup_adult_model = function(model, b, q, s,
                        dispersal_opts = list(),
                        bionomic_opts = list(),
                        eip=15){
  UseMethod("setup_adult_model", model$Mpar)
}
