#' Update States for the `BQS` model: one time step
#'
#' @param t the current time
#' @param model a model defined as a compound [list]
#'
#' @return model, a compound [list]
#' @export
adult_dynamics.BQS = function(t, model){
  with(model,{
    with(c(Mpar, Mvars, terms), {
      eggs_t = ova*psiQ*Q
      # compute variables
      B_t = Mbb %*% B + Mbq %*% Q + Mbs %*% S + Mbl %*% Lambda
      Q_t = Mqb %*% B + Mqq %*% Q + Mqs %*% S
      S_t = Msb %*% B + Msq %*% Q + Mss %*% S + Msl %*% Lambda
      # update variables
      model$Mvars$B = B_t
      model$Mvars$Q = Q_t
      model$Mvars$S = S_t
      model$terms$eggs = eggs_t
    return(model)
  })})
}

#' Store the state variables for the BQS model
#'
#' @param t the current time
#' @param model a model defined as a compound [list]
#'
#' @return a [list] with the most recent states
#' @export
save_states_M.BQS = function(t, model){
  model$states$M$B_t[[t+1]] = model$Mvars$B
  model$states$M$Q_t[[t+1]] = model$Mvars$Q
  model$states$M$S_t[[t+1]] = model$Mvars$S
  model$states$M$eggs_t[[t+1]] = model$Mvars$eggs
  return(model)
}

#' The sum of squared differences between the state and putative steady state
#'
#' @inheritParams compute_diffs_M
#'
#' @return a [numeric] value, the sum of squared differences
#' @export
compute_diffs_M.BQS = function(model){with(model,{
  dfs = sum((Mvars$B - steady$M$B)^2)
  dfs = dfs + sum((Mvars$Q - steady$M$Q)^2)
  dfs = dfs + sum((Mvars$S - steady$M$S)^2)
  return (dfs)
})}

#' Save the state variables in a vector
#'
#' @inheritParams init_states_M
#'
#' @return a [list] with the most recent states
#' @export
init_states_M.BQS = function(model){
  M = list()
  M$B_t    = list()
  M$Q_t    = list()
  M$S_t    = list()
  M$eggs_t = list()
  model$states$M <- M
  model <- save_states_M(0, model)
  return(model)
}

#' Set initial values for the BQS model
#'
#' @inheritParams init_adult_model
#'
#' @return model, a compound [list]
#' @export
init_adult_model.BQS = function(model, M0_opts){
  model = init_adult_model_BQS(model, M0_opts)
  return(model)
}

#' Set initial values for the BQ model
#'
#' @param model a model defined as a compound [list]
#' @param opts a list of options to overwrite defaults
#' @param B0 initial values for blood feeding mosquito density
#' @param Q0 initial values for egg laying mosquito density
#' @param S0 initial values for sugar feeding mosquito density
#'
#' @return model, a compound [list]
#' @export
init_adult_model_BQS = function(model, opts=list(), B0=10, Q0=10, S0=10){with(opts,{
  model$Mvars$B = rep(B0, model$nb)[1:model$nb]
  model$Mvars$Q = rep(Q0, model$nq)[1:model$nq]
  model$Mvars$S = rep(S0, model$ns)[1:model$ns]
  model$Mvars$eggs = rep(0, model$nq)[1:model$nq]
  return(model)
})}


#' Dispersal parameters for the BQS model
#'
#' @param model a model defined as a compound [list]
#' @param opts a [list] of values that overwrite the defaults
#' @param kFb a function that weights points by distance for blood searching
#' @param kFq a function that weights points by distance for habitat searching
#' @param kFs a function that weights points by distance for sugar searching
#' @param wb search weights for blood feeding sites
#' @param wq search weights for habitats
#' @param ws search weights for sugar sites
#' @param stayB the fraction of blood feeding bouts that stay
#' @param stayQ the fraction of egg laying bouts that stay
#' @param stayS the fraction of sugar feeding bouts that stay
#'
#' @return the model, a compound [list]
#' @export
setup_dispersal_BQS = function(model, opts = list(),
                           kFb, kFq, kFs,
                           wb=1, wq=1, ws=1,
                           stayB=0.1, stayQ=0.1, stayS=0.1){
  with(model,{with(opts,{
    model$Mpar$setup$kFb = kFb
    model$Mpar$setup$kFq = kFq
    model$Mpar$setup$kFs = kFs
    model$Mpar$setup$wb = wb
    model$Mpar$setup$wq = wq
    model$Mpar$setup$ws = ws
    model$Mpar$setup$stayB=stayB
    model$Mpar$setup$stayQ=stayQ
    model$Mpar$setup$stayS=stayS
    model = make_Psi_BQS(model)
    return(model)
  })})
}

#' Make dispersal matrices for the BQS model
#'
#' @param model a model defined as a compound [list]
#'
#' @return model, a compound [list]
#' @export
make_Psi_BQS = function(model){
  with(model,{with(Mpar,{with(setup,{
    #from b
    model$Mpar$Psi_bb = make_Psi_xx(b, kFb, wb, stayB)
    model$Mpar$Psi_qb = make_Psi_xy(b, q, kFq, wq)
    model$Mpar$Psi_sb = make_Psi_xy(b, s, kFs, ws)
    #from q
    model$Mpar$Psi_bq = make_Psi_xy(q, b, kFb, wb)
    model$Mpar$Psi_qq = make_Psi_xx(q, kFq, wq, stayQ)
    model$Mpar$Psi_sq = make_Psi_xy(q, s, kFs, ws)
    #from s
    model$Mpar$Psi_bs = make_Psi_xy(s, b, kFb, wb)
    model$Mpar$Psi_qs = make_Psi_xy(s, q, kFq, wq)
    model$Mpar$Psi_ss = make_Psi_xx(s, kFs, ws, stayS)

    return(model)
  })})})
}

#' Bionomic parameters set for the BQ model
#'
#' @param model a model defined as a compound [list]
#' @param opts a [list] of values that overwrite the defaults
#' @param pB the probability of surviving a blood feeding bout
#' @param pQ the probability of surviving an egg laying bout
#' @param pS the probability of surviving a sugar feeding bout
#' @param psiB blood feeding success
#' @param psiQ egg laying success
#' @param psiS sugar feeding success
#' @param sigb transition from blood to sugar feeding
#' @param sigq transition from egg laying to sugar feeding
#' @param sigf transition from egg laying to sugar feeding
#' @param sigL transition from emergence to sugar feeding
#' @param ova female eggs laid during a successful bout, per female
#'
#' @return the parameters, a compound [list]
#' @export
setup_bionomics_BQS = function(model, opts=list(),
                               pB=0.96, pQ=0.96, pS=0.96,
                               psiB=0.9, psiQ=0.9, psiS=0.9,
                               sigb=0.1, sigq=0.1, sigf=0.1,
                               sigL=0.1, ova=20){
  with(opts,{
    model$Mpar$pB = pB
    model$Mpar$pQ = pQ
    model$Mpar$pS = pS
    model$Mpar$psiB = psiB
    model$Mpar$psiQ = psiQ
    model$Mpar$psiS = psiS
    model$Mpar$sigb = sigb
    model$Mpar$sigq = sigq
    model$Mpar$sigf = sigf
    model$Mpar$sigL = sigL
    model$Mpar$ova = ova
  return(model)
})}


#' Make the demographic matrices for the BQS model
#'
#' @param model a model defined as a compound [list]
#'
#' @return model, a compound [list]
#' @export
make_demography_BQS = function(model){
  with(model,{with(Mpar,{with(setup,{
    # from b
    model$Mpar$Mbb = Psi_bb %*% diag(pB*(1-sigb)*(1-psiB), nb)
    model$Mpar$Mqb = Psi_qb %*% diag(pB*psiB, nb)
    model$Mpar$Msb = Psi_sb %*% diag(pB*sigb*psiB, nb)
    # "hardened" adults from q
    model$Mpar$Mbq = Psi_bq %*% diag(pQ*(1-sigf)*psiQ, nq)
    model$Mpar$Mqq = Psi_qq %*% diag(pQ*(1-sigq)*(1-psiQ), nq)
    model$Mpar$Msq = Psi_sq %*% diag(pQ*(sigf*psiQ + sigq*(1-psiQ)), nq)
    # from s
    model$Mpar$Mbs = Psi_bs %*% diag(pS*psiS, ns)
    model$Mpar$Mqs = 0*t(model$Mpar$Msq)
    model$Mpar$Mss = Psi_ss %*% diag(pS*(1-psiS), ns)
    # recently emerged adults
    model$Mpar$Mbl = Psi_bq %*% diag(pQ*(1-sigL), nq)
    model$Mpar$Msl = Psi_sq %*% diag(pQ*sigL, nq)
    # the "hardened adult" mosquito population dispersal matrix
    model$Mpar$bigM = with(model$Mpar, rbind(
      cbind(Mbb, Mbq, Mbs),
      cbind(Mqb, Mqq, Mqs),
      cbind(Msb, Msq, Mss)
    ))
    return(model)
})})})}

#' Setup a BQS model for adult mosquitoes
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
setup_adult_model.BQS = function(model, b, q, s,
                        dispersal_opts = list(),
                        bionomic_opts =list(),
                        eip=15){

  model = setup_dispersal_BQS(model, dispersal_opts)
  model = setup_bionomics_BQS(model, bionomic_opts)
  model = make_demography_BQS(model)
  model$Mpar$eip = eip

  return(model)
}


