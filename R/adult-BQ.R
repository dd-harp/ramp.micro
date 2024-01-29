
#' Simulate adult dynamics for the `BQ` model
#'
#' @inheritParams adult_dynamics
#'
#' @return model, a compound [list]
#' @export
adult_dynamics.BQ = function(t, model){
  with(model,
    with(c(Mpar, Mvars, terms), {
      # compute variables
      eggs_t = ova*psiQ*Q
      B_t = Mbb %*% B + Mbq %*% Q + Mlb %*% Lambda
      Q_t = Mqb %*% B + Mqq %*% Q
      # update variables
      model$Mvars$B = B_t
      model$Mvars$Q = Q_t
      model$terms$eggs = eggs_t
      # return the model
      return(model)
}))}

#' Save the state variables in a vector
#'
#' @param t the current time
#' @param model a model defined as a compound [list]
#'
#' @return a [list] with the most recent states
#' @export
save_states_M.BQ = function(t, model){
  model$states$M$B_t[[t+1]] = model$Mvars$B
  model$states$M$Q_t[[t+1]] = model$Mvars$Q
  model$states$M$eggs_t[[t+1]] = model$Mvars$eggs
  return(model)
}

#' The some of squared differences between two sets of variables
#'
#' @inheritParams compute_diffs_M
#'
#' @return a [numeric] value, the sum of squared differences
#' @export
compute_diffs_M.BQ = function(Mvars1, Mvars2){
  dfs = sum((Mvars1$B - Mvars2$B)^2) + sum((Mvars1$Q - Mvars2$Q)^2)
  return (dfs)
}

#' Save the state variables in a vector
#'
#' @inheritParams init_states_M
#'
#' @return a [list] with the most recent states
#' @export
init_states_M.BQ = function(model){
  M = list()
  M$B_t    = list()
  M$Q_t    = list()
  M$eggs_t = list()
  model$states$M = M
  model = save_states_M(0, model)
  return(model)
}

#' Set initial values for the BQ model
#'
#' @inheritParams init_adult_model
#'
#' @return model, a compound [list]
#' @export
init_adult_model.BQ = function(model, M0_opts){
  model = init_adult_model_BQ(model, M0_opts)
  return(model)
}

#' Set initial values for the BQ model
#'
#' @param model a model defined as a compound [list]
#' @param opts a list of options to overwrite defaults
#' @param B0 initial values for blood feeding mosquito density
#' @param Q0 initial values for egg laying mosquito density
#'
#' @return model, a compound [list]
#' @export
init_adult_model_BQ = function(model, opts=list(), B0=10, Q0=10){with(opts,{
  model$Mvars$B = rep(B0, model$nb)[1:model$nb]
  model$Mvars$Q = rep(Q0, model$nq)[1:model$nq]
  model$Mvars$eggs = rep(0, model$nq)[1:model$nq]
  return(model)
})}


#' Dispersal parameters for the BQ model
#'
#' @param model a model defined as a compound [list]
#' @param opts a [list] of values that overwrite the defaults
#' @param kFb a function that weights points by distance for blood searching
#' @param kFq a function that weights points by distance for habitat searching
#' @param wb search weights for blood feeding sites
#' @param wq search weights for habitats
#' @param stayB the fraction of blood feeding bouts that stay
#' @param stayQ the fraction of egg laying bouts that stay
#'
#' @return the model, a compound [list]
#' @export
setup_dispersal_BQ = function(model, opts = list(),
                              kFb, kFq,
                              wb=1, wq=1,
                              stayB=0.1, stayQ=0.1){
  with(model,{with(opts,{
    model$Mpar$setup = list()
    model$Mpar$setup$kFb = kFb
    model$Mpar$setup$kFq = kFq
    model$Mpar$setup$wb = wb
    model$Mpar$setup$wq = wq
    model$Mpar$setup$stayB=stayB
    model$Mpar$setup$stayQ=stayQ
    model = make_Psi_BQ(model)
    return(model)
  })})
}

#' Make dispersal matrices for the BQS model
#'
#' @param model a model defined as a compound [list]
#'
#' @return model, a compound [list]
#' @export
make_Psi_BQ = function(model){
  with(model,{with(Mpar,{with(setup,{

    #from b
    model$Mpar$Psi_bb = make_Psi_xx(b, kFb, wb, stayB)
    model$Mpar$Psi_qb = make_Psi_xy(b, q, kFq, wq)
    #from q
    model$Mpar$Psi_bq = make_Psi_xy(q, b, kFb, wb)
    model$Mpar$Psi_qq = make_Psi_xx(q, kFq, wq, stayQ)

    return(model)
  })})})
}

#' Bionomic parameters set for the BQ model
#'
#' @param model a model defined as a compound [list]
#' @param opts a [list] of values that overwrite the defaults
#' @param pB the probability of surviving a blood feeding bout
#' @param pQ the probability of surviving an egg laying bout
#' @param psiB blood feeding success
#' @param psiQ egg laying success
#' @param ova female eggs laid during a successful bout, per female
#'
#' @return the parameters, a compound [list]
#' @export
setup_bionomics_BQ = function(model, opts=list(),
                              pB=0.96, pQ=0.96,
                              psiB=0.9, psiQ=0.9,
                              ova=20){
  with(opts,{
    model$Mpar$setup$pB = pB
    model$Mpar$setup$pQ = pQ
    model$Mpar$setup$psiB = psiB
    model$Mpar$psiQ = psiQ
    model$Mpar$ova = ova
    return(model)
  })}


#' Make the demographic matrices for the BQ model
#'
#' @param model a model defined as a compound [list]
#'
#' @return the model
#' @return model, a compound [list]
#' @export
make_demography_BQ = function(model){
  with(model,{with(Mpar,{with(setup,{
    # from b
    model$Mpar$Mbb = Psi_bb %*% diag(pB*(1-psiB), nb)
    model$Mpar$Mqb = Psi_qb %*% diag(pB*psiB, nb)
    # "hardened" adults from q
    model$Mpar$Mbq = Psi_bq %*% diag(pQ*psiQ, nq)
    model$Mpar$Mqq = Psi_qq %*% diag(pQ*(1-psiQ), nq)
    # recently emerged adults
    model$Mpar$Mlb = Psi_bq %*% diag(pQ, nq)
    # the "hardened adult" mosquito population dispersal matrix
    model$Mpar$bigM = with(model$Mpar,
      rbind(
        cbind(Mbb, Mbq),
        cbind(Mqb, Mqq))
      )
    return(model)
})})})}

#' Setup a BQ model for adult mosquitoes
#'
#' @param model a model defined as a compound [list]
#' @param b a point set defining blood feeding sites
#' @param q a point set defining egg laying sites
#' @param s a point set defining sugar feeding sites
#' @param dispersal_opts a [list] to overwrite defaults
#' @param bionomic_opts a [list] to overwrite defaults
#' @param eip the extrinsic incubation period
#'
#' @return a [list] defining a BQ-class adult model
#' @export
setup_adult_model.BQ = function(model, b, q, s=c(),
                           dispersal_opts = list(),
                           bionomic_opts =list(),
                           eip=15){

  model = setup_dispersal_BQ(model, dispersal_opts)
  model = setup_bionomics_BQ(model, bionomic_opts)
  model = make_demography_BQ(model)
  model$Mpar$eip = eip

  return(model)
}
