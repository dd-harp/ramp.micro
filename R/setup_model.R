
#' Setup a behavioral state, micro-simulation model
#'
#' @param b a point set defining blood feeding sites
#' @param q a point set defining egg laying sites
#' @param s a point set defining sugar feeding sites, with a default null value
#' @param kFb a kernel shape for blood searching
#' @param kFq a kernal shape for aquatic habitat searching
#' @param kFs a kernel shape for sugar site searching
#' @param Mname the adult model name
#' @param Lname the aquatic model name
#' @param bionomic_opts a [list] to overwrite defaults
#' @param aquatic_opts a [list] to overwrite defaults
#' @param M0_opts options to overwrite defaults
#' @param L0_opts options to overwrite defaults
#' @param eip the extrinsic incubation period
#'
#' @return a [list] defining a BQ-class adult model
#' @export
setup_model = function(b, q, s=c(), kFb, kFq, kFs = NULL,
                       Mname="BQ", Lname="basicL",
                       bionomic_opts = list(),
                       aquatic_opts = list(),
                       M0_opts = list(),
                       L0_opts = list(),
                       eip=15){

  model = list()

  model$b=b
  model$q=q
  if(!is.null(s)) model$s=s

  model$nb = length(b[,1])
  model$nq = length(q[,1])
  if(!is.null(s)) model$ns = length(s[,1])

  dispersal_opts = list(kFb=kFb, kFq=kFq)
  if(!is.null(kFs)) dispersal_opts$kFs = kFs

  Mpar = list()
  class(Mpar) <- Mname
  model$Mpar = Mpar
  model$Mpar$setup = list()

  model = setup_adult_model(model, b, q, s, dispersal_opts, bionomic_opts, eip)

  model$Mvars = list()
  model = init_adult_model(model, M0_opts)
  model$terms = list()

  Lpar = list()
  class(Lpar) <- Lname
  model$Lpar = Lpar

  model = setup_aquatic_model(model, aquatic_opts)

  model$Lvars = list()
  model = init_aquatic_model(model, L0_opts)

  return(model)
}


#' Basic analysis
#' @param model a model defined as a compound [list]
#' @param burn the burn argument for [steady_state]
#' @param Tx block size argument for [steady_state]
#' @param tol tolerance for [steady_state]
#' @param Tmax runtime for [makeKGV]
#'
#' @return model, a compound [list]
#' @export
basic_analysis = function(model, burn=200, Tx=50, tol =1e-3, Tmax=200){
  model = steady_state(model, burn, Tx, tol)
  model = makeKGV(model, Tmax)
  model = make_tiles(model)
  model = make_all_graphs(model)
  return(model)
}

#' Compute net dispersal matrix to blood feed within one feeding cycle
#'
#' @param model a model defined as a compound [list]
#' @param Tmax the number of time steps
#'
#' @return the model, a compound [list]
#' @export
makeKGV = function(model, Tmax=200){
  model$KGV <- list()
  model = make_Kbq(model, Tmax)
  model = make_Kqb(model, Tmax)
  model = make_Kbb(model)
  model = make_Kqq(model)
  model = compute_G(model, Tmax)
  model = compute_GG(model)
  model = compute_V(model, Tmax)
  model = compute_VC(model)
  return(model)
}

#' Setup a square lattice model with one point set offset and nested within the other
#'
#' @param N the size of the big lattice
#' @param kFb a blood feeding dispersal kernel
#' @param kFq an egg laying dispersal kernel
#' @param q_outside a [logical] switch: if true, habitats define the big lattice
#'
#' @return a model as a compound [list]
#' @export
make_model_squareLattice = function(N, kFb, kFq, q_outside=TRUE){
  x0 = (N-1)/2
  x1 = (N-2)/2
  if(q_outside == TRUE){
    q0 = lattice(N, -x0, x0)
    b0 = lattice(N-1, -x1, x1)
  } else {
    b0 = lattice(N, -x0, x0)
    q0 = lattice(N-1, -x1, x1)
  }
  lattice_mod = setup_model(b0, q0, kFb, kFq, Mname = "BQ")
  lattice_mod = steady_state(lattice_mod)
  lattice_mod = steady_state(lattice_mod)
  lattice_mod$tot = with(lattice_mod$steady_state,c(sum(B), sum(Q)))
  lattice_mod$means = with(lattice_mod$steady_state,c(mean(B), mean(Q)))
  lattice_mod$cv = with(lattice_mod$steady_state,c(cv(B), cv(Q)))
  return(lattice_mod)
}

#' Setup a model on a random uniform set of points
#'
#' @param N the size of the big lattice
#' @param kFb a blood feeding dispersal kernel
#' @param kFq an egg laying dispersal kernel
#' @param q_outside a [logical] switch: if true, habitats define the big lattice
#'
#' @return a model as a compound [list]
#' @export
make_model_unif=function(N, kFb, kFq, q_outside=TRUE){
  x0 = (N-1)/2
  x1 = (N-2)/2
  if(q_outside == TRUE){
    q0 = unif_xy(N^2, -x0, x0)
    b0 = unif_xy((N-1)^2, -x1, x1)
  } else {
    b0 = unif_xy(N^2, -x0, x0)
    q0 = unif_xy((N-1)^2, -x1, x1)
  }
  unif_mod = setup_model(b0, q0, kFb, kFq, Mname = "BQ")
  unif_mod = steady_state(unif_mod)
  unif_mod = steady_state(unif_mod)
  unif_mod$tot = with(unif_mod$steady_state,c(sum(B), sum(Q)))
  unif_mod$means = with(unif_mod$steady_state,c(mean(B), mean(Q)))
  unif_mod$cv = with(unif_mod$steady_state,c(cv(B), cv(Q)))
  return(unif_mod)
}
